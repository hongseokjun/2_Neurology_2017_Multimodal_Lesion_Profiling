function PerpendicularGradientCalculator(output_dir, case_prefix, case_num, left_right, nummesh, intracortical_surface_number, subcortical_surface_number, sampling_space, pve_correction, contrast_str)

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


if(pve_correction)
    pve_postfix = '_RI_corrected';
else
    pve_postfix = '_RI';    
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
    cortical_surf_relative_intensity = zeros(intracortical_surface_number+2, numvert);
    cortical_surf_file = [];
    cortical_surf = [];
    
    for i = 1 : intracortical_surface_number + 2
        %% Step1: Read cortical surfaces of each subject with computed normal vector on each vertex
        if(i == 1)
            basename  = [ case_prefix '_' case_num '_gray_surface_' left_right '_'  num2str(nummesh) postfix_surf ];
        elseif(i == 5)
            basename  = [ case_prefix '_' case_num '_white_surface_' left_right '_'  num2str(nummesh) postfix_surf ];
        else
            basename  = [ case_prefix '_' case_num '_intracortical_surface_' num2str(i-1) '_' left_right '_'  num2str(nummesh) postfix_surf ];
        end

        cortical_surf_file{i} = [ output_dir '/' case_num '/surfaces/'  basename postfix_contrast '.obj' ];
        cortical_surf{i} = SurfStatReadSurf1(cortical_surf_file{i});
        cortical_surf_relative_intensity(i, :) = SurfStatReadData1([ output_dir '/' case_num '/measurement/' basename postfix_contrast pve_postfix '.txt' ]);
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
        if(~exist( [ output_dir '/' case_num '/measurement/' basename postfix_contrast '_pg.txt' ], 'file'))
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
                    surf1 = 4; surf2 = 5;
            end

            temp_surf_1 = cortical_surf{surf1};
            temp_surf_2 = cortical_surf{surf2};
            vwise_distance = sqrt(sum(power((temp_surf_1.coord - temp_surf_2.coord), 2)));

            %% Step3: Sample the relative intensity value at two those sampling points
            %% Step4: Compute perpendicular gradient as following: PG = 1/2*[(sampled_RI_below - sampled_RI_above]/(distance)];
            %% File name: [ case_prefix '_' case_num '_intracortical_surface_' num2str(i-1) '_' left_right '_'  num2str(nummesh) postfix_surf '_pg.txt' ]
            %% variable: cortical_surf_p_gradient
            Sampled_RI_1 = cortical_surf_relative_intensity(surf1, :);
            Sampled_RI_2 = cortical_surf_relative_intensity(surf2, :);
            cortical_surf_p_gradient(i, :) = abs(Sampled_RI_2 - Sampled_RI_1)./vwise_distance;
            SurfStatWriteData([ output_dir '/' case_num '/measurement/' basename postfix_contrast '_pg.txt' ], cortical_surf_p_gradient(i, :), 'a');
        else
            disp([ output_dir '/' case_num '/measurement/' basename postfix_contrast '_pg.txt already exists' ]);
        end
    end
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% subcortical WM surfaces %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subcortical_surf_p_gradient = zeros(subcortical_surface_number, numvert);
    subcortical_surf_relative_intensity = zeros(subcortical_surface_number, numvert);
    subcortical_surf_file = [];
    subcortical_surf = [];

    for i = 1 : subcortical_surface_number
        %% Step1: Read cortical surfaces of each subject with computed normal vector on each vertex
        basename  = [ case_prefix '_' case_num '_white_surface_' num2str(i) '_' left_right '_' num2str(nummesh) postfix_surf ];
        subcortical_surf_file{i} = [ output_dir '/' case_num '/surfaces/'  basename postfix_contrast '.obj' ];
        subcortical_surf{i} = SurfStatReadSurf1(subcortical_surf_file{i});
        subcortical_surf_relative_intensity(i, :) = SurfStatReadData1([ output_dir '/' case_num '/measurement/' basename postfix_contrast '_RI.txt' ]);
    end

    for i = 1 : subcortical_surface_number
        basename  = [ case_prefix '_' case_num '_white_surface_' num2str(i) '_' left_right '_' num2str(nummesh) postfix_surf ];

        %% Step2: Locate two sampling points immediately below and above each vertex
        %% Variable: vwise_distance
        if(~exist( [ output_dir '/' case_num '/measurement/' basename postfix_contrast '_pg.txt' ], 'file'))
            switch i
                case 1 % 1st subcortical WM surface ...
                    surf1 = 1; surf2 = 2;
                case 2 % 2nd subcortical WM surface ...
                    surf1 = 1; surf2 = 3;
                case 3 % 3rd subcortical WM surface ...
                    surf1 = 2; surf2 = 3;
            end

            temp_surf_1 = cortical_surf{surf1};
            temp_surf_2 = cortical_surf{surf2};
            vwise_distance = sqrt(sum(power((temp_surf_1.coord - temp_surf_2.coord), 2)));

            %% Step3: Sample the relative intensity value at two those sampling points
            %% Step4: Compute perpendicular gradient as following: PG = 1/2*[(sampled_RI_below - sampled_RI_above]/(distance)];
            %% File name: [ case_prefix '_' case_num '_intracortical_surface_' num2str(i-1) '_' left_right '_'  num2str(nummesh) postfix_surf '_pg.txt' ]
            %% variable: cortical_surf_p_gradient
            Sampled_RI_1 = subcortical_surf_relative_intensity(surf1, :);
            Sampled_RI_2 = subcortical_surf_relative_intensity(surf2, :);
            subcortical_surf_p_gradient(i, :) = abs(Sampled_RI_2 - Sampled_RI_1)./vwise_distance;
            SurfStatWriteData([ output_dir '/' case_num '/measurement/' basename postfix_contrast '_pg.txt' ], subcortical_surf_p_gradient(i, :), 'a');
        else
            disp([ output_dir '/' case_num '/measurement/' basename postfix_contrast '_pg.txt already exists' ]);
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%
    %% GM-WM blurring %%
    %%%%%%%%%%%%%%%%%%%%
    basename  = [ case_prefix '_' case_num '_white_surface_' left_right '_'  num2str(nummesh) postfix_surf ];
    if(~exist( [ output_dir '/' case_num '/measurement/' basename postfix_contrast '_pg_GM_WM.txt' ], 'file'))                
	    temp_surf_1 = cortical_surf{4};
	    temp_surf_2 = subcortical_surf{1};

	    vwise_distance = sqrt(sum(power((temp_surf_1.coord - temp_surf_2.coord), 2)));

	    ri_mri_name	   = [ output_dir '/' case_num '/' space '/' case_prefix '_' case_num '_ri' postfix_contrast '_gm' postfix_surf '.mnc' ];	    
	    [r, s] = system(['volume_object_evaluate -linear ' ri_mri_name ' ' subcortical_surf_file{1} ' /tmp/tmp_' basename '_' case_prefix '_' case_num postfix_contrast '_RI.txt' ]);      
	    temp_data_1 = cortical_surf_relative_intensity(4, :);
	    temp_data_2 = SurfStatReadData1(['/tmp/tmp_' basename '_' case_prefix '_' case_num postfix_contrast '_RI.txt' ]);	    
	    [r, s] = system(['rm -rf /tmp/tmp_' basename '_' case_prefix '_' case_num postfix_contrast '_RI.txt']);
	    cortical_surf_p_gradient_GM_WM = abs(temp_data_1 - temp_data_2)./vwise_distance;
	    SurfStatWriteData([ output_dir '/' case_num '/measurement/' basename postfix_contrast '_pg_GM_WM.txt' ], cortical_surf_p_gradient_GM_WM, 'a');
    else
    	disp([ output_dir '/' case_num '/measurement/' basename postfix_contrast '_pg_GM_WM.txt already exists' ]);
    end    
end

