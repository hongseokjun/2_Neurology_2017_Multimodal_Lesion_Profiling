function RelativeIntensitySampler(output_dir, case_prefix, case_num, left_right, nummesh, intracortical_surface_number, subcortical_surface_number, sampling_space, contrast_str, overwrite, method)

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
    postfix_contrast = [ '_' contrast_str{j} ];
    postfix_contrast2 = [ '_t1' ];

    if(exist([ output_dir case_num '/measurement/' case_prefix '_' case_num '_white_surface_' num2str(subcortical_surface_number) '_' left_right '_' num2str(nummesh) postfix_surf  postfix_contrast '_RI.txt' ], 'file') ~= 0 && overwrite ~= 1)
	disp(['Already processed: ' case_prefix '_' case_num '_surfaces_' left_right '_' num2str(nummesh) postfix_surf  postfix_contrast '_RI.txt']);
	continue;
    end

    relativeT1_GM_name = [ output_dir case_num '/' space '/' case_prefix '_' case_num '_ri' postfix_contrast '_gm' postfix_surf '.mnc' ];
    relativeT1_WM_name = [ output_dir case_num '/' space '/' case_prefix '_' case_num '_ri' postfix_contrast '_wm' postfix_surf '.mnc' ];
    [r, s] = system(['gunzip -vf ' relativeT1_GM_name '.gz']);
    [r, s] = system(['gunzip -vf ' relativeT1_WM_name '.gz']);
    
    relativeT1_GM = SurfStatReadVol1(relativeT1_GM_name);
    relativeT1_WM = SurfStatReadVol1(relativeT1_WM_name);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% 1) sample the RI values from the volume and map them onto intracortical surfaces
    cortical_surf_relative_intensity = zeros(intracortical_surface_number+2, numvert);
    deformed_cortical_surf_relative_intensity = zeros(intracortical_surface_number, numvert);
    for i = 1 : intracortical_surface_number + 2
        if(i == 1)
            basename  = [ case_prefix '_' case_num '_gray_surface_' left_right '_'  num2str(nummesh) postfix_surf ];
        elseif(i == 5)
            basename  = [ case_prefix '_' case_num '_white_surface_' left_right '_'  num2str(nummesh) postfix_surf  ];
        else
            basename  = [ case_prefix '_' case_num '_intracortical_surface_' num2str(i-1) '_' left_right '_'  num2str(nummesh) postfix_surf ];
        end

        cortical_surf_file{i} = [ output_dir case_num '/surfaces/'  basename  postfix_contrast2 '.obj' ];
	if(method == 1)
	        cortical_surf{i} = SurfStatReadSurf(cortical_surf_file{i});
		cortical_surf_relative_intensity(i, :) = SurfStatVol2Surf(relativeT1_GM, cortical_surf{i});
		SurfStatWriteData( [ output_dir case_num '/measurement/' basename postfix_contrast '_RI.txt' ], cortical_surf_relative_intensity(i, :), 'a');
	elseif(method == 2)
		[r, s] = system(['volume_object_evaluate -linear ' relativeT1_GM_name ' ' cortical_surf_file{i} ' ' output_dir case_num '/measurement/' basename postfix_contrast '_RI.txt']);
	end
    end

    for i = 1 : intracortical_surface_number
        basename  = [ case_prefix '_' case_num '_intracortical_surface_' num2str(i) '_' left_right '_'  num2str(nummesh) '_deformed' postfix_surf ];
        cortical_surf_file_deform{i} = [ output_dir case_num '/surfaces/'  basename postfix_contrast2 '.obj'   ];

        if(exist(cortical_surf_file_deform{i}, 'file') == 2)
            cortical_surf_deform{i} = SurfStatReadSurf(cortical_surf_file_deform{i});
	    if(method == 1)
            	deformed_cortical_surf_relative_intensity(i, :) = SurfStatVol2Surf(relativeT1_GM, cortical_surf_deform{i});
	        SurfStatWriteData( [ output_dir case_num '/measurement/' basename postfix_contrast '_RI.txt' ], deformed_cortical_surf_relative_intensity(i, :), 'a');
	    elseif(method == 2)
	    	[r, s] = system(['volume_object_evaluate -linear ' relativeT1_GM_name ' ' cortical_surf_file_deform{i} ' ' output_dir case_num '/measurement/' basename postfix_contrast '_RI.txt']);
	    end
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% 2) sample the RI values from the volume and map them onto subcortical WM surfaces
    subcortical_surf_relative_intensity = zeros(subcortical_surface_number, numvert);
    for i = 1 : subcortical_surface_number
        basename  = [ case_prefix '_' case_num '_white_surface_' num2str(i) '_' left_right '_' num2str(nummesh) postfix_surf ];
        subcortical_surf_file{i} = [ output_dir case_num '/surfaces/'  basename postfix_contrast2 '.obj' ];
	
	if(method == 1)
        	subcortical_surf{i} = SurfStatReadSurf(subcortical_surf_file{i});
	        subcortical_surf_relative_intensity(i, :) = SurfStatVol2Surf(relativeT1_WM, subcortical_surf{i});
        	SurfStatWriteData( [ output_dir case_num '/measurement/' basename postfix_contrast '_RI.txt' ], subcortical_surf_relative_intensity(i, :), 'a');
	elseif(method == 2)
		[r, s] = system(['volume_object_evaluate -linear ' relativeT1_WM_name ' ' subcortical_surf_file{i} ' ' output_dir case_num '/measurement/' basename postfix_contrast '_RI.txt']);
	end	
    end
        
    [r, s] = system(['gzip -vf ' relativeT1_GM_name ]);
    [r, s] = system(['gzip -vf ' relativeT1_WM_name ]);
end
