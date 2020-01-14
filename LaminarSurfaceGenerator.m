function LaminarSurfaceGenerator(case_prefix, case_num, left_right, intracortical_surface_number, ...
                                 output_dir, mesh_num, roughness_threshold, weight_Jt, weight_Jd, UseLightVersion)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% [INTRODUCTION]
%% This script is to generate intracortical surfaces. 
%% We position vertices of each intra cortical surface at fixed distance 
%% between GM-CSF and GM-WM boundary and then refine the mesh configuration 
%% of the surface by giving uniform intravertex spacing and deformation magnitude constraint. 
%%
%% [IUPUT ARGUMENT]
%% This script presumes that the CIVET processing is already finished for the
%% case. It needs following arguments:
%%      1) case_prefix: e.g. TLE | MCD
%%      2) case_num: e.g 0321_1
%%      3) left_right: left | right
%%      4) intracortical_surface_number: (default := 3)
%%      5) output_dir: The directory which the result will be saved.
%%      6) mesh_num: total number of mesh (81920 | 327680)
%%      7) roughness_threshold: As the first processing step in this pipeline,
%%         roughness of the GM-CSF surface will be measured and the area where
%%         show high roughness will be detected, which is considered as the
%%         spiky shape of surface that need to be locally smoothed before
%%         generating surfaces. (default := 0.2)
%%      8) weight_Jt: A weight coefficient for the intervertex spacing
%%         constraint in energy function. (default := 0.5)
%%      9) weight_Jd: A weight coefficient for the deformation magnitude
%%         constraint in energy function. (default := 0.5)
%%     10) UseLightVersion: If surface deformation is not necessary, then 1, which skip the process of surface deformation
%% 
%% [STEPS]
%%      1) Read surface-related files
%%      2) Calculate mean and standard deviation of intervertex distance in GM-WM and GM-CSF surface
%%      3) Perform Intra-cortical surface positioning
%%      4) Locally deform surfaces, while enforcing a vertex spacing uniformity constraint and deformation magnitude constraint
%%      5) Generate stream lines linking a corresponding vertex between GM-CSF and GM-WM interfaces
%%
%% [OUTPUT]
%% 	measurment:
%%		${PREFIX}_${ID}_${LEFT_RIGHT}_WM_intervertex_dist_avg_${NUMMESH}.txt
%%		${PREFIX}_${ID}_${LEFT_RIGHT}_GM_intervertex_dist_avg_${NUMMESH}.txt
%%		${PREFIX}_${ID}_${LEFT_RIGHT}_[1|2|3]_Jt_${NUMMESH}.txt
%%		${PREFIX}_${ID}_${LEFT_RIGHT}_[1|2|3]_Jd_${NUMMESH}.txt
%%		${PREFIX}_${ID}_${LEFT_RIGHT}_[1|2|3]_intervertex_dist_avg_${NUMMESH}.txt
%%		${PREFIX}_${ID}_${LEFT_RIGHT}_[1|2|3]_intervertex_dist_std_${NUMMESH}.txt
%%		${PREFIX}_${ID}_${LEFT_RIGHT}_[1|2|3]_L_${NUMMESH}.txt
%%
%%	temp:
%%		${PREFIX}_${ID}_${LEFT_RIGHT}_${NUMMESH}_[deformed|]_line.obj
%%		${PREFIX}_${ID}_${LEFT_RIGHT}_${NUMMESH}_[deformed|]_line.vtk	
%%	
%%	surfaces:
%%		${PREFIX}_${ID}_intracortical_surface_[1|2|3]_${LEFT_RIGHT}_${NUMMESH}[_deformed].obj
%%		${PREFIX}_${ID}_intracortical_surface_[1|2|3]_${LEFT_RIGHT}_${NUMMESH}[_deformed].vtk
%%
%% [REFERENCE]
%% 1) For the measurement of surface roughness, please see Lavoue G, 
%%    ACM Trans Applied Perception (2009)
%% 2) For the refinement of intracortical surfaces and energy function
%%    configuration, please see JS Kim, Evans AC, et al NIMG (2005)
%%
%% 
%% Last updated: 12 Mar, 2012 by SeokJun Hong (hong.seok.jun@gmail.com)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1) Reading surface-related files ...
GM_surf_file = [ output_dir '/' case_num '/surfaces/' case_prefix '_' case_num '_gray_surface_' left_right '_' num2str(mesh_num) '_t1.obj'  ];
WM_surf_file = [ output_dir '/' case_num '/surfaces/' case_prefix '_' case_num '_white_surface_' left_right '_' num2str(mesh_num) '_t1.obj' ];

[r, s] = system(['gunzip -vf ' output_dir '/' case_num '/final/' case_prefix '_' case_num '_t1_final.mnc.gz']);
T1  = SurfStatReadVol1([ output_dir '/' case_num '/final/' case_prefix '_' case_num '_t1_final.mnc' ]);
[r, s] = system(['gzip -vf ' output_dir '/' case_num '/final/' case_prefix '_' case_num '_t1_final.mnc']);

WM_surf = SurfStatReadSurf1(WM_surf_file);
GM_surf = SurfStatReadSurf1(GM_surf_file);
WM_surf.edg=SurfStatEdg(WM_surf);
GM_surf.edg=SurfStatEdg(GM_surf);

%% 2) Calculating mean and standard deviation of intervertex distance in GM-WM and GM-CSF surface ...
if(UseLightVersion == 0)
    edg = WM_surf.edg;
    WM_intervertex_dist = zeros(1, size(WM_surf.coord, 2));
    disp('Calculate intervertex distance of GM and WM surface...');
    if(exist([ output_dir '/' case_num '/measurement/' case_prefix '_' case_num '_' left_right '_GM_intervertex_dist_avg_' num2str(mesh_num) '_t1.txt' ], 'file') == 0)
        tstart = tic;
        for i = 1 : size(WM_surf.coord, 2)
            if(mod(i,1000) == 1)
                disp(['WM, percent done: ' num2str(round(i/size(WM_surf.coord, 2)*100),2) '%'])
            end
            NghborEdg=[edg(find(edg(:,1)==i),2)' edg(find(edg(:,2)==i),1)'];
            Dist=distancePoints3d(WM_surf.coord(:,i)', WM_surf.coord(:,NghborEdg)');
            WM_intervertex_dist(i)=mean(Dist);
        end

        edg = GM_surf.edg;
        GM_intervertex_dist = zeros(1, size(GM_surf.coord, 2));
        for i = 1 : size(GM_surf.coord, 2)
            if(mod(i,1000) == 1)
                disp(['GM, percent done: ' num2str(round(i/size(GM_surf.coord, 2)*100),2) '%'])
            end
            NghborEdg=[edg(find(edg(:,1)==i),2)' edg(find(edg(:,2)==i),1)'];
            Dist=distancePoints3d(GM_surf.coord(:,i)', GM_surf.coord(:,NghborEdg)');
            GM_intervertex_dist(i)=mean(Dist);
        end
        telapsed = toc(tstart);
        disp(['Done, elasped time: ' num2str(telapsed/60) ' mins']);

        SurfStatWriteData([ output_dir '/' case_num '/measurement/' case_prefix '_' case_num '_' left_right '_WM_intervertex_dist_avg_' num2str(mesh_num) '_t1.txt' ], WM_intervertex_dist);
        SurfStatWriteData([ output_dir '/' case_num '/measurement/' case_prefix '_' case_num '_' left_right '_GM_intervertex_dist_avg_' num2str(mesh_num) '_t1.txt' ], GM_intervertex_dist);
    else
        WM_intervertex_dist = SurfStatReadData([ output_dir '/' case_num '/measurement/' case_prefix '_' case_num '_' left_right '_WM_intervertex_dist_avg_' num2str(mesh_num) '_t1.txt' ]);
        GM_intervertex_dist = SurfStatReadData([ output_dir '/' case_num '/measurement/' case_prefix '_' case_num '_' left_right '_GM_intervertex_dist_avg_' num2str(mesh_num) '_t1.txt' ]);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%    Intracortical surface vertex generation     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% This algorithm for generating the intermediate surfaces followed the
%% procedure of (JS Kim at al. NIMG 2005) for computing the pial surface
%% reconstruction from the white matter surface. Specifically we start with
%% the white matter surface generated by CIVET and then position the
%% intracortical surface vertices at fixed relative distances along the
%% line segment connecting the gray/white boundary to the pial surface. In
%% addition, we locally deform surfaces, while enforcing a vertex spacing
%% uniformity constraint and deformation magnitude constraint.

%% Term 'Jt': encourage a uniform spacing of vertices 
%% Term 'Jd': give a penalty to an excessive vertex movement from tlink
temp_surf = surfGetNeighbors(GM_surf); 
if(UseLightVersion == 0)
    temp_surf.intervertex_dist = zeros(1, size(GM_surf.coord, 2));
    temp_surf.intervertex_dist_std = zeros(1, size(GM_surf.coord, 2));
    temp_surf.Jt = [];
    temp_surf.Jd = [];
    temp_surf.L =[];
end

intra_cortical_surface(1) = temp_surf;
step_scale_set = [ 0.25 0.125 ];
diff_obj_func  = [ 1000 500 ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Intra-cortical surface positioning
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ------------------------------------ 1-- GM-CSF
%% ------------------------------------ 2-- 1th intra-cortical surface
%% ------------------------------------ 3-- 2nd intra-cortical surface
%% ------------------------------------ 4-- 3rd intra-cortical surface
%% ------------------------------------ 5-- GM-WM

if(exist([ output_dir '/' case_num '/surfaces/' case_prefix '_' case_num '_intracortical_surface_' ...
           num2str(intracortical_surface_number) '_' left_right '_' num2str(mesh_num) '_t1.obj' ], 'file') == 0)
    tstart = tic;
    for i = 1 : intracortical_surface_number
        disp(['Initiation of intracortical surface generation 1: positioning vertices of intracortical surface ' num2str(i) 'th equidistantly']);
        intra_cortical_surface(i) = temp_surf;
        intra_cortical_surface(i).coord = GM_surf.coord - (GM_surf.coord - WM_surf.coord)*(i/(intracortical_surface_number+1)); %% Position vertices of each intra-cortical surface from GM-CSF to GM-WM

        if(UseLightVersion == 0)
            if(exist([ output_dir '/' case_num '/measurement/' case_prefix '_' case_num '_' left_right '_' num2str(intracortical_surface_number) '_L_' num2str(mesh_num) '_t1.txt' ], 'file') == 0)
                disp('Initiation of intracortical surface generation 2: calculating Jn, Jt and intervertex distance...');
                intra_cortical_surface(i).L = mean(sqrt(sum((intra_cortical_surface(i).coord(:,intra_cortical_surface(i).edg(:, 1))' - intra_cortical_surface(i).coord(:,intra_cortical_surface(i).edg(:, 2))').^2, 2)));

                edg = intra_cortical_surface(i).edg;
                for j = 1 : size(GM_surf.coord, 2)
                    if(mod(j,1000) == 1)
                        disp([num2str(i) ' surface, percent done: ' num2str(round(j/size(GM_surf.coord, 2)*100),2) '%'])
                    end

                    NghborEdg=[edg(edg(:,1)==j, 2)' edg(edg(:,2)==j, 1)'];
                    Dist=distancePoints3d(intra_cortical_surface(i).coord(:,j)', intra_cortical_surface(i).coord(:,NghborEdg)');
                    intra_cortical_surface(i).Jt(j) = sum(((Dist - intra_cortical_surface(i).L)/intra_cortical_surface(i).L).^2);

                    %% Point-Line Distance-3Dimensional (see http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html)
                    point_line_distance = sqrt(sum(cross((GM_surf.coord(:, j) - WM_surf.coord(:, j)), (intra_cortical_surface(i).coord(:, j) - WM_surf.coord(:, j))).^2)) / ...
                        sqrt(sum((GM_surf.coord(:, j) - WM_surf.coord(:, j)).^2));
                    intra_cortical_surface(i).Jd(j) = point_line_distance;
                    intra_cortical_surface(i).intervertex_dist(j)=mean(Dist);
                    intra_cortical_surface(i).intervertex_dist_std(j)=std(Dist);
                end
                SurfStatWriteData([ output_dir '/' case_num '/measurement/' case_prefix '_' case_num '_' left_right '_' num2str(i) '_Jt_' num2str(mesh_num) '_t1.txt' ], ...
                    intra_cortical_surface(i).Jt);
                SurfStatWriteData([ output_dir '/' case_num '/measurement/' case_prefix '_' case_num '_' left_right '_' num2str(i) '_Jd_' num2str(mesh_num) '_t1.txt' ], ...
                    intra_cortical_surface(i).Jd);
                SurfStatWriteData([ output_dir '/' case_num '/measurement/' case_prefix '_' case_num '_' left_right '_' num2str(i) '_intervertex_dist_avg_' num2str(mesh_num) '_t1.txt' ], ...
                    intra_cortical_surface(i).intervertex_dist);
                SurfStatWriteData([ output_dir '/' case_num '/measurement/' case_prefix '_' case_num '_' left_right '_' num2str(i) '_intervertex_dist_std_' num2str(mesh_num) '_t1.txt' ], ...
                    intra_cortical_surface(i).intervertex_dist_std);
                SurfStatWriteData([ output_dir '/' case_num '/measurement/' case_prefix '_' case_num '_' left_right '_' num2str(i) '_L_' num2str(mesh_num) '_t1.txt' ], ...
                    intra_cortical_surface(i).L);
            else
                intra_cortical_surface(i).Jt = SurfStatReadData([ output_dir '/' case_num '/measurement/' case_prefix '_' case_num '_' left_right '_' num2str(i) '_Jt_' num2str(mesh_num) '_t1.txt' ]);
                intra_cortical_surface(i).Jd = SurfStatReadData([ output_dir '/' case_num '/measurement/' case_prefix '_' case_num '_' left_right '_' num2str(i) '_Jd_' num2str(mesh_num) '_t1.txt' ]);
                intra_cortical_surface(i).intervertex_dist = SurfStatReadData([ output_dir '/' case_num '/measurement/' case_prefix '_' case_num '_' left_right '_' num2str(i) '_intervertex_dist_avg_' num2str(mesh_num) '_t1.txt' ]);
                intra_cortical_surface(i).intervertex_dist_std  = SurfStatReadData([ output_dir '/' case_num '/measurement/' case_prefix '_' case_num '_' left_right '_' num2str(i) '_intervertex_dist_std_' num2str(mesh_num) '_t1.txt' ]);
                intra_cortical_surface(i).L = SurfStatReadData([ output_dir '/' case_num '/measurement/' case_prefix '_' case_num '_' left_right '_' num2str(i) '_L_' num2str(mesh_num) '_t1.txt' ]);
            end
        end
    end
    telapsed = toc(tstart);

    disp(['Initiation done, elasped time: ' num2str(telapsed/60) ' mins']);

    intra_cortical_surface_curr = intra_cortical_surface;

    if(UseLightVersion == 0)
        fid = fopen([ output_dir '/' case_num '/logs/' case_prefix '_' case_num '_' left_right '_object_function_value_' num2str(mesh_num) '_t1.txt'], 'w');
        disp(['Intracortical surface generation start!']);

        tstart = tic;
        for i = 1 : intracortical_surface_number
            energy_func_org = sum(weight_Jt*intra_cortical_surface_curr(i).Jt + weight_Jd*intra_cortical_surface_curr(i).Jd);
            for s = 1 : size(step_scale_set, 2)
                count = 0;
                step_scale = step_scale_set(s);
                edg = intra_cortical_surface(i).edg;
                temp_deformed_surf = intra_cortical_surface_curr(i);
                energy_func_past = 100000000000; energy_func_curr = sum(weight_Jt*temp_deformed_surf.Jt + weight_Jd*temp_deformed_surf.Jd);

                while(abs(energy_func_past - energy_func_curr) > diff_obj_func(s))
                    count = count + 1;
                    energy_func_past = energy_func_curr;
                    for j = 1 : size(GM_surf.coord, 2)
                        temp_deformed_surf_for_this_iter = temp_deformed_surf;
                        NghborEdg=[edg(edg(:,1)==j, 2)' edg(edg(:,2)==j, 1)'];
                        triangle_index = temp_deformed_surf_for_this_iter.tri((temp_deformed_surf_for_this_iter.tri(:, 1) == j), :);
                        triangle_index = [ triangle_index; temp_deformed_surf_for_this_iter.tri((temp_deformed_surf_for_this_iter.tri(:, 2) == j), :) ];
                        triangle_index = [ triangle_index; temp_deformed_surf_for_this_iter.tri((temp_deformed_surf_for_this_iter.tri(:, 3) == j), :) ];

                        Jt_curr = [ temp_deformed_surf_for_this_iter.Jt(j) temp_deformed_surf_for_this_iter.Jt(NghborEdg) ];
                        Jd_curr = [ temp_deformed_surf_for_this_iter.Jd(j) temp_deformed_surf_for_this_iter.Jd(NghborEdg) ];
                        J_curr = sum(weight_Jt*Jt_curr + weight_Jd*Jd_curr);

                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %% deformed vector calculation ...
                        %% It will generate (2 * the number of nearest neighbours)
                        %% vectors ...
                        deformed_vector = [];

                        for k = 1 : size(triangle_index, 1)
                            two_neighbours = triangle_index(k, (triangle_index(k, :)~=j));
                            deformed_vector_temp = ((temp_deformed_surf_for_this_iter.coord(:, two_neighbours(1)) - temp_deformed_surf_for_this_iter.coord(:, j)) + ...
                                (temp_deformed_surf_for_this_iter.coord(:, two_neighbours(2)) - temp_deformed_surf_for_this_iter.coord(:, j))) / 2;

                            normal_deformed_vector =  dot(deformed_vector_temp, double(temp_deformed_surf_for_this_iter.normal(:, j)))*(temp_deformed_surf_for_this_iter.normal(:, j));
                            tangential_deformed_vector = deformed_vector_temp - normal_deformed_vector;
                            deformed_vector = [ deformed_vector tangential_deformed_vector ];
                        end
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                        Jt_temp = zeros(size(NghborEdg, 2)+1, size(deformed_vector, 2));
                        Jd_temp = zeros(size(NghborEdg, 2)+1, size(deformed_vector, 2));
                        for k = 1 : size(deformed_vector, 2)
                            temp_deformed_surf_for_this_iter.coord(:, j) = temp_deformed_surf.coord(:, j) + step_scale*deformed_vector(:, k);
                            point_line_distance = sqrt(sum(cross((GM_surf.coord(:, j) - WM_surf.coord(:, j)), (temp_deformed_surf_for_this_iter.coord(:, j) - WM_surf.coord(:, j))).^2)) / ...
                                sqrt(sum((GM_surf.coord(:, j) - WM_surf.coord(:, j)).^2));
                            intravertex_dist_limit = abs(0.5*((i/(intracortical_surface_number+1))*(GM_intervertex_dist(j)-WM_intervertex_dist(j))+ WM_intervertex_dist(j)));

                            if(intravertex_dist_limit > point_line_distance)
                                %% re-calculation of 'Jn', 'Jt' and intervertex distance terms ...
                                Dist=kron(ones(1, size(temp_deformed_surf_for_this_iter.coord(:,NghborEdg), 2)), temp_deformed_surf_for_this_iter.coord(:,j)) - temp_deformed_surf_for_this_iter.coord(:,NghborEdg);
                                Dist = sqrt(sum(Dist.^2));
                                Jt_temp(1, k) = sum(((Dist - temp_deformed_surf_for_this_iter.L)/temp_deformed_surf_for_this_iter.L).^2);
                                Jd_temp(1, k) = point_line_distance;

                                for n = 1 : size(NghborEdg, 2)
                                    NghborEdg_temp=[edg(edg(:,1)==NghborEdg(n), 2)' edg(edg(:,2)==NghborEdg(n), 1)'];

                                    %% re-calculation of 'Jn', 'Jt' and intervertex distance terms of neighbours ...
                                    Dist = kron(ones(1, size(temp_deformed_surf_for_this_iter.coord(:,NghborEdg_temp), 2)), temp_deformed_surf_for_this_iter.coord(:,NghborEdg(n))) - ...
                                        temp_deformed_surf_for_this_iter.coord(:,NghborEdg_temp);
                                    Dist = sqrt(sum(Dist.^2));
                                    Jt_temp(1+n, k) = sum(((Dist - temp_deformed_surf_for_this_iter.L)/temp_deformed_surf_for_this_iter.L).^2);
                                    Jd_temp(1+n, k) = temp_deformed_surf_for_this_iter.Jd(NghborEdg(n));
                                end
                            else
                                Jt_temp(:, k) = 1000;
                                Jd_temp(:, k) = 1000;
                            end
                        end

                        J_temp = sum(weight_Jt*Jt_temp + weight_Jd*Jd_temp, 1);
                        [J_min_temp, J_min_temp_index] = min(J_temp);

                        if(J_min_temp >= J_curr)
                            temp_deformed_surf_for_this_iter.coord(:, j) = temp_deformed_surf.coord(:, j);
                        elseif(J_min_temp < J_curr)
                            temp_deformed_surf_for_this_iter.coord(:, j) = temp_deformed_surf.coord(:, j) + step_scale*deformed_vector(:, J_min_temp_index);

                            %% new normal vector after vertex movement ...
                            u1=temp_deformed_surf_for_this_iter.coord(:,triangle_index(:,1));
                            d1=temp_deformed_surf_for_this_iter.coord(:,triangle_index(:,2))-u1;
                            d2=temp_deformed_surf_for_this_iter.coord(:,triangle_index(:,3))-u1;
                            c=cross(d1,d2,1);
                            temp_deformed_surf_for_this_iter.normal(:,j) = mean(c, 2)/sqrt(sum(mean(c, 2).^2));

                            %% re-calculation of 'Jn', 'Jt' and intervertex distance terms ...
                            Dist=kron(ones(1, size(temp_deformed_surf_for_this_iter.coord(:,NghborEdg), 2)), temp_deformed_surf_for_this_iter.coord(:,j)) - temp_deformed_surf_for_this_iter.coord(:,NghborEdg);
                            Dist = sqrt(sum(Dist.^2));
                            temp_deformed_surf_for_this_iter.Jt(j) = sum(((Dist - temp_deformed_surf_for_this_iter.L)/temp_deformed_surf_for_this_iter.L).^2);
                            point_line_distance = sqrt(sum(cross((GM_surf.coord(:, j) - WM_surf.coord(:, j)), (temp_deformed_surf_for_this_iter.coord(:, j) - WM_surf.coord(:, j))).^2)) / ...
                                sqrt(sum((GM_surf.coord(:, j) - WM_surf.coord(:, j)).^2));
                            temp_deformed_surf_for_this_iter.Jd(j) = point_line_distance;

                            for n = 1 : size(NghborEdg, 2)
                                NghborEdg_temp=[edg(edg(:,1)==NghborEdg(n), 2)' edg(edg(:,2)==NghborEdg(n), 1)'];
                                triangle_index_temp = temp_deformed_surf_for_this_iter.tri((temp_deformed_surf_for_this_iter.tri(:, 1) == NghborEdg(n)), :);
                                triangle_index_temp = [ triangle_index_temp; temp_deformed_surf_for_this_iter.tri((temp_deformed_surf_for_this_iter.tri(:, 2) == NghborEdg(n)), :) ];
                                triangle_index_temp = [ triangle_index_temp; temp_deformed_surf_for_this_iter.tri((temp_deformed_surf_for_this_iter.tri(:, 3) == NghborEdg(n)), :) ];

                                %% new normal vector of neighbours after vertex movement ...
                                u1=temp_deformed_surf_for_this_iter.coord(:,triangle_index_temp(:,1));
                                d1=temp_deformed_surf_for_this_iter.coord(:,triangle_index_temp(:,2))-u1;
                                d2=temp_deformed_surf_for_this_iter.coord(:,triangle_index_temp(:,3))-u1;
                                c=cross(d1,d2,1);
                                temp_deformed_surf_for_this_iter.normal(:, NghborEdg(n)) = mean(c, 2)/sqrt(sum(mean(c, 2).^2));

                                %% re-calculation of 'Jn', 'Jt' and intervertex distance terms of neighbours ...
                                Dist = kron(ones(1, size(temp_deformed_surf_for_this_iter.coord(:,NghborEdg_temp), 2)), temp_deformed_surf_for_this_iter.coord(:,NghborEdg(n))) - ...
                                    temp_deformed_surf_for_this_iter.coord(:,NghborEdg_temp);
                                Dist = sqrt(sum(Dist.^2));
                                temp_deformed_surf_for_this_iter.Jt(NghborEdg(n)) = sum(((Dist - temp_deformed_surf_for_this_iter.L)/temp_deformed_surf_for_this_iter.L).^2);
                            end
                            energy_func_curr = sum(weight_Jt*temp_deformed_surf_for_this_iter.Jt + weight_Jd*temp_deformed_surf_for_this_iter.Jd);
                            disp(['iteration: ' num2str(count) ' surface ' num2str(i) ' step scale: ' num2str(step_scale) ', vertex ' num2str(j) 'th energy function : ' num2str(energy_func_curr)]);
                        end
                        temp_deformed_surf = temp_deformed_surf_for_this_iter;
                        energy_func_curr = sum(weight_Jt*temp_deformed_surf.Jt + weight_Jd*temp_deformed_surf.Jd);
                    end

                    temp_deformed_surf.L = mean(sqrt(sum((temp_deformed_surf.coord(:,temp_deformed_surf.edg(:, 1))' - temp_deformed_surf.coord(:,temp_deformed_surf.edg(:, 2))').^2, 2)));
                    u1=temp_deformed_surf.coord(:,temp_deformed_surf.tri(:,1));
                    d1=temp_deformed_surf.coord(:,temp_deformed_surf.tri(:,2))-u1;
                    d2=temp_deformed_surf.coord(:,temp_deformed_surf.tri(:,3))-u1;
                    c=cross(d1,d2,1);
                    temp_deformed_surf.normal=zeros(3,size(temp_deformed_surf.coord, 2));
                    for j=1:3
                        for k=1:3
                            temp_deformed_surf.normal(k,:)=temp_deformed_surf.normal(k,:)+accumarray(temp_deformed_surf.tri(:,j),c(k,:)')';
                        end
                    end
                    temp_deformed_surf.normal=temp_deformed_surf.normal./(ones(3,1)*sqrt(sum(temp_deformed_surf.normal.^2,1)));
                    fprintf(fid, 'Intracortical layer: %d, step scale: %f, iteration: %d, difference of energy function between original and deformed: %f --> %f\n', ...
                        i, step_scale, count, energy_func_org, energy_func_curr);
                end
                intra_cortical_surface_curr(i) = temp_deformed_surf;
            end
            SurfStatWriteData([ output_dir '/' case_num '/measurement/' case_prefix '_' case_num '_' left_right '_' num2str(i) '_Jt_' num2str(mesh_num) '_t1.txt' ], ...
                intra_cortical_surface_curr(i).Jt);
            SurfStatWriteData([ output_dir '/' case_num '/measurement/' case_prefix '_' case_num '_' left_right '_' num2str(i) '_Jd_' num2str(mesh_num) '_t1.txt' ], ...
                intra_cortical_surface_curr(i).Jd);
            SurfStatWriteData([ output_dir '/' case_num '/measurement/' case_prefix '_' case_num '_' left_right '_' num2str(i) '_intervertex_dist_avg_' num2str(mesh_num) '_t1.txt' ], ...
                intra_cortical_surface_curr(i).intervertex_dist);
            SurfStatWriteData([ output_dir '/' case_num '/measurement/' case_prefix '_' case_num '_' left_right '_' num2str(i) '_intervertex_dist_std_' num2str(mesh_num) '_t1.txt' ], ...
                intra_cortical_surface_curr(i).intervertex_dist_std);
            SurfStatWriteData([ output_dir '/' case_num '/measurement/' case_prefix '_' case_num '_' left_right '_' num2str(i) '_L_' num2str(mesh_num) '_t1.txt' ], ...
                intra_cortical_surface_curr(i).L);
        end
        telapsed = toc(tstart);
        disp(['Intracortical surface generation done, elasped time: ' num2str(telapsed/60) ' mins']);
        fclose(fid);
    end
else
    for i = 1 : intracortical_surface_number
        temp = SurfStatReadSurf1([ output_dir '/' case_num '/surfaces/' case_prefix '_' case_num '_intracortical_surface_' num2str(i) '_' left_right '_' num2str(mesh_num) '_t1.obj' ]);
        temp.edg = temp_surf.edg;
        temp.nbr = temp_surf.nbr;
        if(UseLightVersion == 0)
            temp.intervertex_dist = SurfStatReadData([ output_dir '/' case_num '/measurement/' case_prefix '_' case_num '_' left_right '_' num2str(i) '_intervertex_dist_avg_' num2str(mesh_num) '_t1.txt' ]);
            temp.intervertex_dist_std  = SurfStatReadData([ output_dir '/' case_num '/measurement/' case_prefix '_' case_num '_' left_right '_' num2str(i) '_intervertex_dist_std_' num2str(mesh_num) '_t1.txt' ]);
            temp.Jt = SurfStatReadData([ output_dir '/' case_num '/measurement/' case_prefix '_' case_num '_' left_right '_' num2str(i) '_Jt_' num2str(mesh_num) '_t1.txt' ]);
            temp.Jd = SurfStatReadData([ output_dir '/' case_num '/measurement/' case_prefix '_' case_num '_' left_right '_' num2str(i) '_Jd_' num2str(mesh_num) '_t1.txt' ]);
            temp.L = SurfStatReadData([ output_dir '/' case_num '/measurement/' case_prefix '_' case_num '_' left_right '_' num2str(i) '_L_' num2str(mesh_num) '_t1.txt' ]);
        end
        intra_cortical_surface_curr(i) = temp;
        intra_cortical_surface(i) = temp;
    end
end

intra_cortical_surface_deformed = intra_cortical_surface_curr;

%% 3) Generate stream lines linking a corresponding vertex between GM-CSF and GM-WM interfaces
if(UseLightVersion)
    postfix='';
else
    postfix='_deformed';    
end

if(exist([ output_dir '/' case_num '/temp/' case_prefix '_' case_num '_' left_right '_' num2str(mesh_num) postfix '_stream_line.vtk' ], 'file') == 0)
    streams = cell(1, size(GM_surf.coord, 2));
    for i = 1 : size(GM_surf.coord, 2)
        streams{i} = [ WM_surf.coord(:, i)' ];
        for j = 1 : intracortical_surface_number
            streams{i} = [ streams{i}; intra_cortical_surface_deformed(j).coord(:, i)' ];
        end
        streams{i} = [ streams{i}; GM_surf.coord(:, i)' ];
    end

    st_lines = ExploreDTI_to_tractStructure(streams,T1.dim(1:3)',T1.vox');
    for f = 1 : st_lines.nFiberNr
        st_lines.fiber(f).data.cellIndex_mean = f;
        goodcoords = find(~isnan(st_lines.fiber(f).xyzFiberCoord(:,1)));

        st_lines.fiber(f).xyzFiberCoord = st_lines.fiber(f).xyzFiberCoord(goodcoords,:);
        st_lines.fiber(f).nFiberLength = length(goodcoords);
        st_lines.fiber(f).nSelectFiberEndPoint = length(goodcoords)-1;
    end
    sampleLines = tracts_removeData(st_lines);
    status = f_save_obj_lines(sampleLines, [ output_dir '/' case_num '/temp/' case_prefix '_' case_num '_' left_right '_' num2str(mesh_num) postfix '_stream_line.obj' ]);
    save_tract_vtk(sampleLines, [ output_dir '/' case_num '/temp/' case_prefix '_' case_num '_' left_right '_' num2str(mesh_num) postfix '_stream_line.vtk' ], 'ASCII');
end
    
for i = 1 : intracortical_surface_number
    if(exist([ output_dir '/' case_num '/surfaces/' case_prefix '_' case_num '_intracortical_surface_' ...
               num2str(intracortical_surface_number) '_' left_right '_' num2str(mesh_num) '.obj' ], 'file') == 0)
           SurfStatWriteSurf1([ output_dir '/' case_num '/surfaces/' case_prefix '_' case_num '_intracortical_surface_' ...
                                num2str(i) '_' left_right '_' num2str(mesh_num) '_t1.obj' ], intra_cortical_surface(i));
    end
                   
    if(UseLightVersion==0)                   
        SurfStatWriteSurf1([ output_dir '/' case_num '/surfaces/' case_prefix '_' case_num '_intracortical_surface_' ...
                             num2str(i) '_' left_right '_' num2str(mesh_num) '_deformed_t1.obj' ], intra_cortical_surface_deformed(i));   
    end
end
