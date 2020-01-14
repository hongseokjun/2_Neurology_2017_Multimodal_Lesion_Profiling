function PerpendicularGradientCalculator_absolute_distance(process_dir, outdir, case_prefix, case_num, left_right, nummesh, intracortical_surface_number, subcortical_surface_number, sampling_space, contrast_str)

if(nummesh == 81920)
    numvert = 40962;
elseif(nummesh == 327680)
    numvert = 163842;
end

if(strcmp(sampling_space, 'native'))
    postfix_surf='_native';
    space = 'native';
elseif(strcmp(sampling_space, 'tal'))
    postfix_surf='';
    space = 'final';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Perpendicular gradient analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1) Locate two sampling points geometrically immediately below and above each vertex
%% 2) Sample relative intensity value at those two sampling points
%% 3) Compute perpendicular gradient using the following equation: PG = [(sampled_RI_below - sampled_RI_above]/(distance of two sample points)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% intra-cortical surfaces %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j = 1 : size(contrast_str, 2)
    
    postfix_contrast = [ '_' contrast_str{j} ];
    
    cortical_surf_p_gradient = zeros(intracortical_surface_number+2, numvert);
    cortical_surf_file = [];
    cortical_surf = [];
    
    relativeT1_GM_name = [ process_dir case_num '/' space '/' case_prefix '_' case_num '_ri' postfix_contrast '_gm' postfix_surf '.mnc' ];
    relativeT1_WM_name = [ process_dir case_num '/' space '/' case_prefix '_' case_num '_ri' postfix_contrast '_wm' postfix_surf '.mnc' ];
    [r, s] = system(['gunzip -vf ' relativeT1_GM_name '.gz']);
    [r, s] = system(['gunzip -vf ' relativeT1_WM_name '.gz']);
    
    relativeT1_GM = SurfStatReadVol1(relativeT1_GM_name);
    relativeT1_WM = SurfStatReadVol1(relativeT1_WM_name);
    
    for i = 1 : intracortical_surface_number + 2 + 1
        %% Step1: Read cortical surfaces of each subject with computed normal vector on each vertex
        if(i == 1)
            basename  = [ case_prefix '_' case_num '_gray_surface_' left_right '_'  num2str(nummesh) postfix_surf ];
        elseif(i == 5)
            basename  = [ case_prefix '_' case_num '_white_surface_' left_right '_'  num2str(nummesh) postfix_surf ];
        elseif(i == 6)
            basename  = [ case_prefix '_' case_num '_white_surface_1_' left_right '_' num2str(nummesh) postfix_surf ];            
        else
            basename  = [ case_prefix '_' case_num '_intracortical_surface_' num2str(i-1) '_' left_right '_'  num2str(nummesh) postfix_surf ];
        end

        cortical_surf_file{i} = [ process_dir '/' case_num '/surfaces/'  basename postfix_contrast '.obj' ];
        cortical_surf{i} = SurfStatReadSurf1(cortical_surf_file{i});
%         cortical_surf_relative_intensity(i, :) = SurfStatReadData1([ process_dir '/' case_num '/measurement/' basename postfix_contrast pve_postfix '.txt' ]);
    end

    for i = 1 : intracortical_surface_number + 2
        %% Step1: Read cortical surfaces of each subject with computed normal vector on each vertex
        if(i == 1)
            basename  = [ case_prefix '_' case_num '_gray_surface_' left_right '_'  num2str(nummesh) postfix_surf ];
        elseif(i == 5)
            basename  = [ case_prefix '_' case_num '_white_surface_' left_right '_'  num2str(nummesh) postfix_surf ];
        else
            basename  = [ case_prefix '_' case_num '_intracortical_surface_' num2str(i-1) '_' left_right '_'  num2str(nummesh) postfix_surf ];
        end

        %% Step2: Locate two sampling points immediately below and above each vertex
        %% Variable: vwise_distance
        if(~exist( [ outdir '/' basename postfix_contrast '_pg_abs_dist.txt' ], 'file'))
            switch i
                case 1 % GM-CSF boundary ...
                    surf1 = 1; surf2 = 2;
                case 2 % Intracortical surface 1 ...
                    surf1 = 1; surf2 = 3;
                case 3 % Intracortical surface 2 ...
                    surf1 = 2; surf2 = 4;
                case 4 % Intracortical surface 3 ...
                    surf1 = 3; surf2 = 5;
                case 5 % WM-GM boundary ...
                    surf1 = 4; surf2 = 6;
            end

            temp_surf_1 = cortical_surf{surf1};
            temp_surf_2 = cortical_surf{surf2};
            curr_surf   = cortical_surf{i};
            diffvec_1 = temp_surf_1.coord-curr_surf.coord;
            diffvec_2 = temp_surf_2.coord-curr_surf.coord;            
            norm_diffvec_1 = diffvec_1./repmat(sqrt(sum(diffvec_1.^2, 1)), 3, 1);
            norm_diffvec_2 = diffvec_2./repmat(sqrt(sum(diffvec_2.^2, 1)), 3, 1);
            
            above_curr_surf = curr_surf;
            below_curr_surf = curr_surf;
            
            above_curr_surf.coord = curr_surf.coord + norm_diffvec_1*0.5;
            below_curr_surf.coord = curr_surf.coord + norm_diffvec_2*0.5;
            
            SurfStatWriteSurf1([outdir '/' basename postfix_contrast '_above.obj'], above_curr_surf);
            SurfStatWriteSurf1([outdir '/' basename postfix_contrast '_below.obj'], below_curr_surf);
            
            [r, s] = system(['volume_object_evaluate -linear ' relativeT1_GM_name ' ' [outdir '/' basename postfix_contrast '_above.obj'] ' ' outdir '/' basename postfix_contrast '_RI_above.txt']);
            [r, s] = system(['volume_object_evaluate -linear ' relativeT1_GM_name ' ' [outdir '/' basename postfix_contrast '_below.obj'] ' ' outdir '/' basename postfix_contrast '_RI_below.txt']);
            
            Sampled_RI_1 = SurfStatReadData1([outdir '/' basename postfix_contrast '_RI_above.txt']);
            Sampled_RI_2 = SurfStatReadData1([outdir '/' basename postfix_contrast '_RI_below.txt']);
            Sampled_RI_1(isnan(Sampled_RI_1)) = 0;
            Sampled_RI_2(isnan(Sampled_RI_2)) = 0;
            
            %% Step3: Sample the relative intensity value at two those sampling points
            %% Step4: Compute perpendicular gradient as following: PG = 1/2*[(sampled_RI_below - sampled_RI_above]/(distance)];
            %% File name: [ case_prefix '_' case_num '_intracortical_surface_' num2str(i-1) '_' left_right '_'  num2str(nummesh) postfix_surf '_pg.txt' ]
            %% variable: cortical_surf_p_gradient            
            cortical_surf_p_gradient(i, :) = abs(Sampled_RI_2 - Sampled_RI_1);
            SurfStatWriteData([ outdir '/' basename postfix_contrast '_pg_abs_dist.txt' ], cortical_surf_p_gradient(i, :), 'a');
        else
            disp([ outdir '/' basename postfix_contrast '_pg.txt already exists' ]);
        end        
    end
    delete([outdir '/*.obj']);
    delete([outdir '/*_above.txt']);
    delete([outdir '/*_below.txt']);
end

