function RelativeIntensitySampler(output_dir, dti_dir, case_prefix, case_num, left_right, nummesh, intracortical_surface_number, subcortical_surface_number, sampling_space, contrast_str, overwrite, method)

if(nummesh == 81920)
    numvert = 40962;
elseif(nummesh == 327680)
    numvert = 163842;
end

if(strcmp(sampling_space, 'native'))
    space='native';
    postfix_surf='_native';
elseif(strcmp(sampling_space, 'tal'))
    space='final';
    postfix_surf='';
end

for j = 1 : size(contrast_str, 2)
    postfix_contrast = [ '_' contrast_str{j} '_bbr' ];

    if(exist([ output_dir case_num '/measurement/' case_prefix '_' case_num '_white_surface_5_' left_right '_' num2str(nummesh) postfix_surf  postfix_contrast '.txt' ], 'file') ~= 0 && overwrite ~= 1)
	disp(['Already processed: ' case_prefix '_' case_num '_surfaces_' left_right '_' num2str(nummesh) postfix_surf  postfix_contrast '.txt']);
	continue;
    end

    DTI_FA_name = [ dti_dir '/' case_num '/mncdti/' case_prefix '_' case_num '_dtifsl_FA.mnc' ];
    DTI_MD_name = [ dti_dir '/' case_num '/mncdti/' case_prefix '_' case_num '_dtifsl_MD.mnc' ];
    [r, s] = system(['gunzip -vf ' DTI_FA_name '.gz']);
    [r, s] = system(['gunzip -vf ' DTI_MD_name '.gz']);
    
    DTI_FA = SurfStatReadVol1(DTI_FA_name);
    DTI_MD = SurfStatReadVol1(DTI_MD_name);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% 1) sample the RI values from the volume and map them onto intracortical surfaces
    cortical_surf_FA = zeros(intracortical_surface_number+2, numvert);
    cortical_surf_MD = zeros(intracortical_surface_number+2, numvert);    
    for i = 1 : intracortical_surface_number + 2
        if(i == 1)
            basename  = [ case_prefix '_' case_num '_gray_surface_' left_right '_'  num2str(nummesh) postfix_surf ];
        elseif(i == 5)
            basename  = [ case_prefix '_' case_num '_white_surface_' left_right '_'  num2str(nummesh) postfix_surf  ];
        else
            basename  = [ case_prefix '_' case_num '_intracortical_surface_' num2str(i-1) '_' left_right '_'  num2str(nummesh) postfix_surf ];
        end

        cortical_surf_file{i} = [ output_dir case_num '/surfaces/'  basename  postfix_contrast '.obj' ];
	if(method == 1)
	        cortical_surf{i} = SurfStatReadSurf(cortical_surf_file{i});
		cortical_surf_FA(i, :) = SurfStatVol2Surf(DTI_FA, cortical_surf{i});
		cortical_surf_MD(i, :) = SurfStatVol2Surf(DTI_MD, cortical_surf{i});
		SurfStatWriteData( [ output_dir case_num '/measurement/' basename postfix_contrast '_FA.txt' ], cortical_surf_FA(i, :), 'a');		
		SurfStatWriteData( [ output_dir case_num '/measurement/' basename postfix_contrast '_MD.txt' ], cortical_surf_MD(i, :), 'a');
	elseif(method == 2)
		[r, s] = system(['volume_object_evaluate -linear ' DTI_FA_name ' ' cortical_surf_file{i} ' ' output_dir case_num '/measurement/' basename postfix_contrast '_FA.txt']);
		[r, s] = system(['volume_object_evaluate -linear ' DTI_MD_name ' ' cortical_surf_file{i} ' ' output_dir case_num '/measurement/' basename postfix_contrast '_MD.txt']);
	end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% 2) sample the RI values from the volume and map them onto subcortical WM surfaces
    subcortical_surf_FA = zeros(subcortical_surface_number, numvert);
    subcortical_surf_MD = zeros(subcortical_surface_number, numvert);
    for i = 1 : subcortical_surface_number
        basename  = [ case_prefix '_' case_num '_white_surface_' num2str(i) '_' left_right '_' num2str(nummesh) postfix_surf ];
        subcortical_surf_file{i} = [ output_dir case_num '/surfaces/'  basename postfix_contrast '.obj' ];
	
	if(method == 1)
        	subcortical_surf{i} = SurfStatReadSurf(subcortical_surf_file{i});
	        subcortical_surf_FA(i, :) = SurfStatVol2Surf(DTI_FA, subcortical_surf{i});
		subcortical_surf_MD(i, :) = SurfStatVol2Surf(DTI_MD, subcortical_surf{i});
        	SurfStatWriteData( [ output_dir case_num '/measurement/' basename postfix_contrast '_FA.txt' ], subcortical_surf_FA(i, :), 'a');
		SurfStatWriteData( [ output_dir case_num '/measurement/' basename postfix_contrast '_MD.txt' ], subcortical_surf_MD(i, :), 'a');
	elseif(method == 2)
		[r, s] = system(['volume_object_evaluate -linear ' DTI_FA_name ' ' subcortical_surf_file{i} ' ' output_dir case_num '/measurement/' basename postfix_contrast '_FA.txt']);
		[r, s] = system(['volume_object_evaluate -linear ' DTI_MD_name ' ' subcortical_surf_file{i} ' ' output_dir case_num '/measurement/' basename postfix_contrast '_MD.txt']);
	end	
    end
    [r, s] = system(['gzip -vf ' DTI_FA_name ]);
    [r, s] = system(['gzip -vf ' DTI_MD_name ]);   
end
