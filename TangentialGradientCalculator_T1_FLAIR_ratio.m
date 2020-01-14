function TangentialGradientCalculator(output_dir, case_prefix, case_num, left_right, nummesh, intracortical_surface_number, subcortical_surface_number, sampling_space, pve_correction, contrast_str)

if(nummesh == 81920)
    numvert = 40962;
elseif(nummesh == 327680)
    numvert = 163842;
end

if(strcmp(sampling_space, 'native'))
    postfix_surf='_native';
elseif(strcmp(sampling_space, 'tal'))
    postfix_surf='';    
end

if(pve_correction)
    pve_postfix = '_RI_corrected';
else
    pve_postfix = '_RI';    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Tangential gradient analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1) Sample RI values from 5 or 6 neighbouring vertices of a given vertex
%% 2) Compute tangential graident using the following equation:
%%    TG(x) = (?(abs(RI(i)-RI(x))/d))/# of neighbourings

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% intra-cortical surfaces %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for c = 1 : size(contrast_str, 2)   
    postfix_contrast = [ '_' contrast_str{c} ];
    postfix_contrast2 = [ '_t1' ];
    
    cortical_surf_t_gradient = zeros(intracortical_surface_number+2, numvert);
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

        cortical_surf_file{i} = [ output_dir '/' case_num '/surfaces/'  basename postfix_contrast2 '.obj' ];
        temp = SurfStatReadSurf1(cortical_surf_file{i});
        cortical_surf{i} = surfGetNeighborsHong(temp);
        cortical_surf_relative_intensity(i, :) = SurfStatReadData1([ output_dir '/' case_num '/measurement/' basename postfix_contrast pve_postfix '.txt' ]);

        if(~exist( [ output_dir '/' case_num '/measurement/' basename postfix_contrast '_tg.txt' ], 'file'))
            %% Step2: Sample RI values from 5 or 6 neighbouring vertices of a given vertex
            %% Step3: Compute tangential graident using the following equation: TG(x) = (?(abs(RI(i)-RI(x))/d))/# of neighbourings
            for j = 1 : numvert
                if(mod(j,1000) == 1)
                    disp([num2str(i) 'th intracortical surface, percent done: ' num2str((j/numvert)*100) '%'])
                end
                neighbourings = cortical_surf{i}.nbr(j, :);
                neighbourings = neighbourings(neighbourings~=0);
                d = sqrt(sum(power(kron(cortical_surf{i}.coord(:, j), ones(1, size(neighbourings, 2)))-cortical_surf{i}.coord(:, neighbourings), 2), 1));
                cortical_surf_t_gradient(i, j) = mean(abs(cortical_surf_relative_intensity(i, neighbourings) - cortical_surf_relative_intensity(i, j))./d);
            end

            %% File name: [ case_prefix '_' case_num '_intracortical_surface_' num2str(i-1) '_' left_right '_'  num2str(nummesh) postfix_surf '_tg.txt' ]
            %% variable: cortical_surf_t_gradient
            SurfStatWriteData([ output_dir '/' case_num '/measurement/' basename postfix_contrast '_tg.txt' ], cortical_surf_t_gradient(i, :), 'a');
        else
            disp([ output_dir '/' case_num '/measurement/' basename postfix_contrast '_tg.txt already exists' ]);
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% subcortical WM surfaces %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subcortical_surf_t_gradient = zeros(subcortical_surface_number, numvert);
    subcortical_surf_relative_intensity = zeros(subcortical_surface_number, numvert);
    subcortical_surf_file = [];
    subcortical_surf = [];

    for i = 1 : subcortical_surface_number
        %% Step1: Read cortical surfaces of each subject with computed normal vector on each vertex
        basename  = [ case_prefix '_' case_num '_white_surface_' num2str(i) '_' left_right '_' num2str(nummesh) postfix_surf ];
        subcortical_surf_file{i} = [ output_dir '/' case_num '/surfaces/'  basename postfix_contrast2 '.obj' ];
        temp = SurfStatReadSurf1(subcortical_surf_file{i});
        subcortical_surf{i} = surfGetNeighborsHong(temp);
        subcortical_surf_relative_intensity(i, :) = SurfStatReadData1([ output_dir '/' case_num '/measurement/' basename postfix_contrast '_RI.txt' ]);

        if(~exist( [ output_dir '/' case_num '/measurement/' basename postfix_contrast '_tg.txt' ], 'file'))
            %% Step2: Sample RI values from 5 or 6 neighbouring vertices of a given vertex
            %% Step3: Compute tangential graident using the following equation: TG(x) = (?(abs(RI(i)-RI(x))/d))/# of neighbourings
            for j = 1 : numvert
                if(mod(j,1000) == 1)
                    disp([num2str(i) 'th subcortical surface, percent done: ' num2str((j/numvert)*100) '%'])
                end
                neighbourings = subcortical_surf{i}.nbr(j, :);
                neighbourings = neighbourings(neighbourings~=0);
                d = sqrt(sum(power(kron(subcortical_surf{i}.coord(:, j), ones(1, size(neighbourings, 2)))-subcortical_surf{i}.coord(:, neighbourings), 2), 1));
                subcortical_surf_t_gradient(i, j) = mean(abs(subcortical_surf_relative_intensity(i, neighbourings) - subcortical_surf_relative_intensity(i, j))./d);
            end

            %% File name: [ case_prefix '_' case_num '_intracortical_surface_' num2str(i-1) '_' left_right '_'  num2str(nummesh) postfix_surf '_tg.txt' ]
            %% variable: cortical_surf_t_gradient
            SurfStatWriteData([ output_dir '/' case_num '/measurement/' basename postfix_contrast '_tg.txt' ], subcortical_surf_t_gradient(i, :), 'a');
        else
            disp([ output_dir '/' case_num '/measurement/' basename postfix_contrast '_tg.txt already exists' ]);
        end
    end
end
