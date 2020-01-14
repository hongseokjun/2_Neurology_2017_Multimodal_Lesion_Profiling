clear all
close all
addpath('/data/noel/noel2/local/matlab/surfstat_chicago');

%% Read pre-calculated zscore of T1 and FLAIR in the lesion
for read_matfiles = 1
    
    patient_comp = 0;
    control_comp = 0;
    disease_control_comp = 1;
    Kernel = 2;
    
    if patient_comp == 1
        load(['zscore_database_perilesional_fMRI_sm_' num2str(Kernel) '.mat']);
        load('individual_surface_multimodal.mat'); % geodesic distance computed based on original surface
                                                   % and modified using small sphere kernel to solve jammed WM surfaces
    end
    
    if control_comp == 1
        load(['zscore_database_perilesional_fMRI_sm_' num2str(Kernel) '_cont.mat']);    
        load('individual_surface_multimodal_cont.mat'); % geodesic distance computed based on original surface
                                                        % and modified using small sphere kernel to solve jammed WM surfaces
    end
    
    if disease_control_comp == 1
        
        load(['zscore_database_perilesional_fMRI_sm_' num2str(Kernel) '_disease_cont.mat']);
        load('individual_surface_multimodal_disease_cont.mat'); % geodesic distance computed based on original surface
                                                                % and modified using small sphere kernel to solve jammed WM surfaces
        
    end    
    
    load([ 'p_value_set_perilesion_profile_T1_FLAIR.mat' ]);
    load([ 'p_value_set_perilesion_profile_DTI.mat' ]);
    load([ 'p_value_set_perilesion_profile_all_modalities.mat' ]);
    
    load('/local_raid/seokjun/01_project/04_IntracorticalAnalysis/01_analysis/pipeline/postprocessing_zscore_visualization/colormap_noel1.mat');
    
end

%% Read demograpic data and set up some parameters and variables
for read_demodata_setup_params = 1
    
    OUTPATH = '/local_raid/seokjun/01_project/04_IntracorticalAnalysis/03_result/lesion_profile/';
    Group_cont                = 'control'
    Prefix_cont               = 'TLE'
    Cases_pat                 = '/local_raid/seokjun/01_project/04_IntracorticalAnalysis/01_analysis/Demographic_data_FCD_fMRI.txt'
    Cases_pat_temp            = '/local_raid/seokjun/01_project/04_IntracorticalAnalysis/01_analysis/Demographic_data_FCD.txt'
    Group_pat                 = 'FCD'
    Prefix_pat                = 'mcd'
    Cases_cont                = '/local_raid/seokjun/01_project/04_IntracorticalAnalysis/01_analysis/Demographic_data_control_fMRI.txt'
    Cases_cont_temp           = '/local_raid/seokjun/01_project/04_IntracorticalAnalysis/01_analysis/Demographic_data_control.txt'
    Cases_disease_cont        = '/local_raid/seokjun/01_project/04_IntracorticalAnalysis/01_analysis/Demographic_data_TLE_final.txt'
    Group_disease_cont        = 'TLE'
    Prefix_disease_cont       = 'TLE'    
    
    Left_Right                = 'both'
    NumIntSurf                = 3
    NumSubSurf                = 3
    NumMesh                   = 81920
    SamplingSpace             = 'native'
    img_contrast              = {'t1'}
    Parametric                = 'quadratic'
    average_surface_dir       = '/local_raid/seokjun/01_project/04_IntracorticalAnalysis/03_result/average_surfaces/'
        
    visualization = 1;
    save_file = 0;
    perilesional_distance = 20; % 20 mm
    pl_feature_read_again = 0;  % 0 = read presaved files | 1 = compute features again

    
    fid = fopen(Cases_cont);
    demo = textscan(fid, '%s%d%s%s%s', 'Delimiter', ',', 'CollectOutput', 1);
    case_num_cont = demo{1};
    age_cont = demo{2};
    gender_cont = demo{3}(:, 1);    
    fclose(fid);
    
    fid = fopen(Cases_cont_temp);
    demo = textscan(fid, '%s%d%s%s%s', 'Delimiter', ',', 'CollectOutput', 1);
    case_num_cont_temp = demo{1};   
    fclose(fid);
        
    case_num_cont_temp(strcmp(case_num_cont_temp, '306_1')) = [];
    case_num_cont_temp(strcmp(case_num_cont_temp, '313_1')) = [];
    case_num_cont_temp(strcmp(case_num_cont_temp, '322_1')) = [];
    [a b] = intersect(case_num_cont_temp, case_num_cont);
    outlier_vector = ones(size(case_num_cont_temp, 1), 1);
    outlier_vector(b) = 0;

    
    if control_comp == 1
        
        geoDist_t1_cont_new(logical(outlier_vector)) = [];
        geoDist_flair_cont_new(logical(outlier_vector)) = [];
        geoDist_dti_cont_new(logical(outlier_vector)) = [];        
        geoDist_fMRI_cont_new(logical(outlier_vector)) = [];        
        
    end
    
    outlier_cases = { '306_1', '313_1', '322_2', '342_1', '343_1' };
    case_cont_real = 1:size(case_num_cont, 1);
    outlier_vector = zeros(size(case_num_cont, 1), 1);
    for i = 1 : size(outlier_cases, 2)
        outlier_vector = outlier_vector + strcmp(case_num_cont, outlier_cases{i});
    end
    case_num_cont(logical(outlier_vector)) = [];
    
    if control_comp == 1
        fMRI_ALFF_z_cont(logical(outlier_vector), :) = [];
        fMRI_ReHo_z_cont(logical(outlier_vector), :) = [];
    end
    
    fid = fopen(Cases_pat);
    demo = textscan(fid, '%s%d%s%s%s', 'Delimiter', ',', 'CollectOutput', 1);
    case_num_pat = demo{1};
    age_pat      = demo{2};
    gender_pat   = demo{3}(:, 1);
    histo_pat    = demo{3}(:, 3);
    
    fid = fopen(Cases_pat_temp);
    demo = textscan(fid, '%s%f%s%d%s%d%d%s%s', 'Delimiter', ',', 'CollectOutput', 1);
    case_num_pat_temp = demo{1}; case_num_pat_temp{strcmp(case_num_pat_temp, '080_1')} = '080_2';
    fclose(fid);
    
    fid = fopen(Cases_disease_cont);
    demo = textscan(fid, '%s%f%s%s', 'Delimiter', ',', 'CollectOutput', 1);
    case_num_disease_cont       = demo{1};
    age_disease_cont            = demo{2};
    gender_disease_cont         = demo{3}(:, 1);
    seizure_lateralization  = demo{3}(:, 2);
    fclose(fid);       
    
    lesion_label_dir = '/local_raid/seokjun/CTFCD-1.2.0_64/Lesion/lesion_surf/';
    lesion_label_data = zeros(length(case_num_pat), 81924);
    for i = 1 : length(case_num_pat)
        lesion_label_data(i, :) = SurfStatReadData( { [ lesion_label_dir '/' Prefix_pat '_' case_num_pat{i} '_label_union_left_rsl.txt' ], ...
            [ lesion_label_dir '/' Prefix_pat '_' case_num_pat{i} '_label_union_right_rsl.txt' ] } );
    end
    
    % remove 061 075 079 from geoDist_dti
    [C, ia, ib] = intersect(case_num_pat_temp, case_num_pat);
    
    if patient_comp == 1
        
        geoDist_fMRI     = geoDist_fMRI(ia);
        geoDist_fMRI_new = geoDist_fMRI_new(ia);
        
    end
    
    surf_template = SurfStatReadSurf({'surf_reg_model_left.obj', 'surf_reg_model_right.obj'});
    mask_template = SurfStatMaskCut(surf_template);
    
%     FDR_threshold_final_perilesion_profiling = FDR([p_value_perilesion_T1_FLAIR p_value_perilesion_DTI], 0.05); %% 0.0025 
    FDR_threshold_final_perilesion_profiling = FDR(p_value_set_all_perilesion_profiling, 0.05); %% 0.0049 
    
end

%% Compute the peri-lesional features along the cortical mantle
for feature_extraction = 1
    
    %% mean_z_lesion
    %% 1   = ALFF of fMRI              -> 1   p1
    %% 2   = ReHo of fMRI              -> 2   p2    
    
    if patient_comp == 1
        
        feature_table = [];
        feature_parameters = [
            { 'fMRI_ALFF_z', 1, 1 }; %% 01 ALFF 1
            { 'fMRI_ReHo_z', 1, 1 }; %% 02 ReHO 2
            ];
        
        if(pl_feature_read_again == 0)
            load(['pl_feature_set_sm_' num2str(Kernel) '_fMRI.mat']);
        elseif(pl_feature_read_again == 1)
            pl_feature_set = [];
            for k = 1 : size(feature_parameters, 1)
                pl_feature_set(k).name = feature_parameters{k, 1};
                pl_feature_set(k).data = zeros(perilesional_distance+1, size(case_num_pat, 1), feature_parameters{k, 2});
            end
            
            lesion_label_dir = '/data/noel/noel6/CTFCD-1.2.0_64/Lesion/lesion_surf/';
            for j = 1 : size(case_num_pat, 1)
                disp(['case ' case_num_pat{j} ' start!']);
                idx = find(strcmp(case_num_pat, case_num_pat{j}));
                if(exist([ lesion_label_dir '/' Prefix_pat '_' case_num_pat{j} '_label_union_left_rsl.txt' ], 'file'))
                    lesion_label_data = SurfStatReadData( { [ lesion_label_dir '/' Prefix_pat '_' case_num_pat{j} '_label_union_left_rsl.txt' ], ...
                        [ lesion_label_dir '/' Prefix_pat '_' case_num_pat{j} '_label_union_right_rsl.txt' ] } );
                    for k = 1 : size(feature_parameters, 1)
                        surf_num = feature_parameters{k, 2};
                        if(strfind(feature_parameters{k, 1}, 'T1'))
                            geoDist = geoDist_t1;
                            if(strfind(feature_parameters{k, 1}, 'Sub'))
                                geoDist = geoDist_t1_new;
                            end
                            disp(feature_parameters{k, 1});
                        elseif(strfind(feature_parameters{k, 1}, 'FLAIR'))
                            geoDist = geoDist_flair;
                            if(strfind(feature_parameters{k, 1}, 'Sub'))
                                geoDist = geoDist_flair_new;
                            end
                            disp(feature_parameters{k, 1});
                        elseif(strfind(feature_parameters{k, 1}, 'fMRI'))
                            geoDist = geoDist_fMRI;
                        end
                        
                        for s = 1 : surf_num
                            if(surf_num > 1)
                                eval([ 'feature = squeeze(' feature_parameters{k, 1} '(s, :, j));' ]);
                            else
                                eval([ 'feature = ' feature_parameters{k, 1} '(j, :);' ]);
                            end
                            
                            %% features
                            pl_feature = mean(feature(logical(lesion_label_data)));
                            surf_start = feature_parameters{k, 3}-1;
                            small_num_vert = [];
                            for pl = 1 : 1 : perilesional_distance
                                vert_idx = find(geoDist{j}.surf(surf_start+s,:) > pl-1 & geoDist{j}.surf(surf_start+s,:) <= pl);
                                if(length(vert_idx) < 5)
                                    small_num_vert = [ small_num_vert pl+1 ];
                                end
                                exclude = find(isnan(feature(vert_idx)) | isinf(feature(vert_idx)) | mask_template(vert_idx) == 0);
                                vert_idx(exclude) = [];
                                pl_feature = [ pl_feature; mean(feature(vert_idx)) ];
                            end
                            
                            %% Remove nan or outlier values
                            %% 1) NAN values can happen when there is no vertex within a given distance range (e.g. 1mm-2mm)
                            %% 2) Outlier can happen when only small number vertices are selected
                            %%    and these vertex fall into the area of noisy values (>3SD or <-3SD)
                            start_nan = 0;
                            end_nan = 0;
                            for pl = 1 : perilesional_distance + 1
                                if(isnan(pl_feature(pl)) & start_nan == 0)
                                    start_nan = pl;
                                elseif(isnan(pl_feature(pl)) & start_nan ~= 0)
                                    end_nan = pl;
                                end
                            end
                            
                            if(end_nan == 0)
                                end_nan = start_nan;
                            end
                            
                            if(start_nan == 0 && end_nan == 0)
                                for sn = 1 : length(small_num_vert)
                                    pl_feature(small_num_vert(sn)) = (pl_feature(small_num_vert(sn)+1)+pl_feature(small_num_vert(sn)-1))/2;
                                end
                            else
                                num_nan = end_nan - start_nan + 1;
                                unit = (pl_feature(end_nan + 1) - pl_feature(start_nan - 1))/(num_nan + 1);
                                for pl = 1 : num_nan
                                    pl_feature(start_nan - 1 + pl) = pl_feature(start_nan-1) + unit;
                                end
                                for sn = 1 : length(small_num_vert)
                                    pl_feature(small_num_vert(sn)) = (pl_feature(small_num_vert(sn)+1)+pl_feature(small_num_vert(sn)-1))/2;
                                end
                            end
                            outlier_idx = find(pl_feature > mean(pl_feature)+3*std(pl_feature) | pl_feature < mean(pl_feature)-3*std(pl_feature));
                            for sn = 1 : length(outlier_idx)
                                if(outlier_idx(sn)==1)
                                    continue;
                                elseif(outlier_idx(sn)==perilesional_distance+1)
                                    pl_feature(outlier_idx(sn)) = pl_feature(outlier_idx(sn)-1);
                                else
                                    pl_feature(outlier_idx(sn)) = (pl_feature(outlier_idx(sn)+1)+pl_feature(outlier_idx(sn)-1))/2;
                                end
                            end
                            pl_feature_set(k).data(:, j, s) = pl_feature;
                        end
                    end
                end
            end
        end
        
        perilesional_distance_org = perilesional_distance;
        pl_feature_set_org = pl_feature_set;
        
    end
    
    if control_comp == 1
        
        feature_table = [];
        feature_parameters = [
            { 'fMRI_ALFF_z_cont', 1, 3 }; %% 01 ALFF 1
            { 'fMRI_ReHo_z_cont', 1, 3 }; %% 02 ReHO 2
            ];
        
        if(pl_feature_read_again == 0)
            load(['pl_feature_set_sm_' num2str(Kernel) '_fMRI_cont.mat']);
        elseif(pl_feature_read_again == 1)
            pl_feature_set = [];
            for k = 1 : size(feature_parameters, 1)
                pl_feature_set(k).name = feature_parameters{k, 1};
                pl_feature_set(k).data = zeros(perilesional_distance+1, size(case_num_cont, 1), feature_parameters{k, 2}, size(case_num_pat, 1));
            end
            
            lesion_label_dir = '/data/noel/noel6/CTFCD-1.2.0_64/Lesion/lesion_surf/';
            for c = 1 : size(case_num_cont, 1)
                
                for k = 1 : size(feature_parameters, 1)
                    surf_num = feature_parameters{k, 2};
                    geoDist = geoDist_fMRI_cont_new;
                    
                    disp(feature_parameters{k, 1});
                    
                    for s = 1 : surf_num
                        if(surf_num > 1)
                            eval([ 'feature = squeeze(' feature_parameters{k, 1} '(s, :, c));' ]);
                        else
                            eval([ 'feature = ' feature_parameters{k, 1} '(c, :);' ]);
                        end
                        
                        for j = 1 : size(case_num_pat_temp, 1)
                            disp(['case :' case_num_cont{c} ', feature : ', feature_parameters{k, 1}, ', surf :', num2str(s) ', pat : ', case_num_pat_temp{j}]);
                            idx = find(strcmp(case_num_pat_temp, case_num_pat_temp{j}));
                            
                            if(exist([ lesion_label_dir '/' Prefix_pat '_' case_num_pat_temp{j} '_label_union_left_rsl.txt' ], 'file'))
                                lesion_label_data = SurfStatReadData( { [ lesion_label_dir '/' Prefix_pat '_' case_num_pat_temp{j} '_label_union_left_rsl.txt' ], ...
                                    [ lesion_label_dir '/' Prefix_pat '_' case_num_pat_temp{j} '_label_union_right_rsl.txt' ] } );
                                
                                %% features
                                pl_feature = mean(feature(logical(lesion_label_data)));
                                surf_start = feature_parameters{k, 3}-1;
                                small_num_vert = [];
                                
                                for pl = 1 : 1 : perilesional_distance
                                    vert_idx = find(geoDist{c}.surf(j, :) > pl-1 & geoDist{c}.surf(j, :) <= pl);
                                    if(length(vert_idx) < 5)
                                        small_num_vert = [ small_num_vert pl+1 ];
                                    end
                                    exclude = find(isnan(feature(vert_idx)) | isinf(feature(vert_idx)) | mask_template(vert_idx) == 0);
                                    vert_idx(exclude) = [];
                                    pl_feature = [ pl_feature; mean(feature(vert_idx)) ];
                                end
                                
                                %% Remove nan or outlier values
                                %% 1) NAN values can happen when there is no vertex within a given distance range (e.g. 1mm-2mm)
                                %% 2) Outlier can happen when only small number vertices are selected
                                %%    and these vertex fall into the area of noisy values (>3SD or <-3SD)
                                start_nan = 0;
                                end_nan = 0;
                                for pl = 1 : perilesional_distance + 1
                                    if(isnan(pl_feature(pl)) & start_nan == 0)
                                        start_nan = pl;
                                    elseif(isnan(pl_feature(pl)) & start_nan ~= 0)
                                        end_nan = pl;
                                    end
                                end
                                
                                if(end_nan == 0)
                                    end_nan = start_nan;
                                end
                                
                                if(start_nan == 0 && end_nan == 0)
                                    for sn = 1 : length(small_num_vert)
                                        pl_feature(small_num_vert(sn)) = (pl_feature(small_num_vert(sn)+1)+pl_feature(small_num_vert(sn)-1))/2;
                                    end
                                else
                                    num_nan = end_nan - start_nan + 1;
                                    unit = (pl_feature(end_nan + 1) - pl_feature(start_nan - 1))/(num_nan + 1);
                                    for pl = 1 : num_nan
                                        pl_feature(start_nan - 1 + pl) = pl_feature(start_nan-1) + unit;
                                    end
                                    for sn = 1 : length(small_num_vert)
                                        pl_feature(small_num_vert(sn)) = (pl_feature(small_num_vert(sn)+1)+pl_feature(small_num_vert(sn)-1))/2;
                                    end
                                end
                                outlier_idx = find(pl_feature > mean(pl_feature)+3*std(pl_feature) | pl_feature < mean(pl_feature)-3*std(pl_feature));
                                for sn = 1 : length(outlier_idx)
                                    if(outlier_idx(sn)==1)
                                        continue;
                                    elseif(outlier_idx(sn)==perilesional_distance+1)
                                        pl_feature(outlier_idx(sn)) = pl_feature(outlier_idx(sn)-1);
                                    else
                                        pl_feature(outlier_idx(sn)) = (pl_feature(outlier_idx(sn)+1)+pl_feature(outlier_idx(sn)-1))/2;
                                    end
                                end
                                pl_feature_set(k).data(:, c, s, j) = pl_feature;
                                
                            end
                        end
                    end
                end 
                
            end
        end
        
        perilesional_distance_org = perilesional_distance;
        pl_feature_set_org = pl_feature_set;
        
    end
    
    if disease_control_comp == 1
        
        pl_feature_read_again = 1;
        feature_table = [];
        feature_parameters = [
            { 'fMRI_ALFF_z_disease_cont', 1, 3 }; %% 01 ALFF 1
            { 'fMRI_ReHo_z_disease_cont', 1, 3 }; %% 02 ReHO 2
            ];
        
        if(pl_feature_read_again == 0)
            load(['pl_feature_set_sm_' num2str(Kernel) '_fMRI_disease_cont.mat']);
        elseif(pl_feature_read_again == 1)
            pl_feature_set = [];
            for k = 1 : size(feature_parameters, 1)
                pl_feature_set(k).name = feature_parameters{k, 1};
                pl_feature_set(k).data = zeros(perilesional_distance+1, size(case_num_disease_cont, 1), feature_parameters{k, 2}, size(case_num_pat, 1));
            end
                        
            for c = 1 : size(case_num_disease_cont, 1)
                
                for k = 1 : size(feature_parameters, 1)
                    surf_num = feature_parameters{k, 2};
                    geoDist = geoDist_fMRI_disease_cont_new;
                    
                    disp(feature_parameters{k, 1});
                    
                    for s = 1 : surf_num
                        if(surf_num > 1)
                            eval([ 'feature = squeeze(' feature_parameters{k, 1} '(s, :, c));' ]);
                        else
                            eval([ 'feature = ' feature_parameters{k, 1} '(c, :);' ]);
                        end
                        
                        for j = 1 : size(case_num_pat_temp, 1)
                            disp(['case :' case_num_disease_cont{c} ', feature : ', feature_parameters{k, 1}, ', surf :', num2str(s) ', pat : ', case_num_pat_temp{j}]);
                            idx = find(strcmp(case_num_pat_temp, case_num_pat_temp{j}));
                            
                            if(exist([ lesion_label_dir '/' Prefix_pat '_' case_num_pat_temp{j} '_label_union_left_rsl.txt' ], 'file'))
                                lesion_label_data = SurfStatReadData( { [ lesion_label_dir '/' Prefix_pat '_' case_num_pat_temp{j} '_label_union_left_rsl.txt' ], ...
                                    [ lesion_label_dir '/' Prefix_pat '_' case_num_pat_temp{j} '_label_union_right_rsl.txt' ] } );
                                
                                %% features
                                pl_feature = mean(feature(logical(lesion_label_data)));
                                surf_start = feature_parameters{k, 3}-1;
                                small_num_vert = [];
                                
                                for pl = 1 : 1 : perilesional_distance
                                    vert_idx = find(geoDist{c}.surf(j, :) > pl-1 & geoDist{c}.surf(j, :) <= pl);
                                    if(length(vert_idx) < 5)
                                        small_num_vert = [ small_num_vert pl+1 ];
                                    end
                                    exclude = find(isnan(feature(vert_idx)) | isinf(feature(vert_idx)) | mask_template(vert_idx) == 0);
                                    vert_idx(exclude) = [];
                                    pl_feature = [ pl_feature; mean(feature(vert_idx)) ];
                                end
                                
                                %% Remove nan or outlier values
                                %% 1) NAN values can happen when there is no vertex within a given distance range (e.g. 1mm-2mm)
                                %% 2) Outlier can happen when only small number vertices are selected
                                %%    and these vertex fall into the area of noisy values (>3SD or <-3SD)
                                start_nan = 0;
                                end_nan = 0;
                                for pl = 1 : perilesional_distance + 1
                                    if(isnan(pl_feature(pl)) & start_nan == 0)
                                        start_nan = pl;
                                    elseif(isnan(pl_feature(pl)) & start_nan ~= 0)
                                        end_nan = pl;
                                    end
                                end
                                
                                if(end_nan == 0)
                                    end_nan = start_nan;
                                end
                                
                                if(start_nan == 0 && end_nan == 0)
                                    for sn = 1 : length(small_num_vert)
                                        pl_feature(small_num_vert(sn)) = (pl_feature(small_num_vert(sn)+1)+pl_feature(small_num_vert(sn)-1))/2;
                                    end
                                else
                                    num_nan = end_nan - start_nan + 1;
                                    unit = (pl_feature(end_nan + 1) - pl_feature(start_nan - 1))/(num_nan + 1);
                                    for pl = 1 : num_nan
                                        pl_feature(start_nan - 1 + pl) = pl_feature(start_nan-1) + unit;
                                    end
                                    for sn = 1 : length(small_num_vert)
                                        pl_feature(small_num_vert(sn)) = (pl_feature(small_num_vert(sn)+1)+pl_feature(small_num_vert(sn)-1))/2;
                                    end
                                end
                                outlier_idx = find(pl_feature > mean(pl_feature)+3*std(pl_feature) | pl_feature < mean(pl_feature)-3*std(pl_feature));
                                for sn = 1 : length(outlier_idx)
                                    if(outlier_idx(sn)==1)
                                        continue;
                                    elseif(outlier_idx(sn)==perilesional_distance+1)
                                        pl_feature(outlier_idx(sn)) = pl_feature(outlier_idx(sn)-1);
                                    else
                                        pl_feature(outlier_idx(sn)) = (pl_feature(outlier_idx(sn)+1)+pl_feature(outlier_idx(sn)-1))/2;
                                    end
                                end
                                pl_feature_set(k).data(:, c, s, j) = pl_feature;
                                
                            end
                        end
                    end
                end 
                
            end
        end
        
        perilesional_distance_org = perilesional_distance;
        pl_feature_set_org = pl_feature_set;
        
    end
    
end

%% Start perilesional analyses
for perilesional_analyses = 1
    
    if patient_comp == 1
        
        for i = 1 : size(pl_feature_set_org, 2)
            pl_feature_set(i).data = pl_feature_set_org(i).data(1:perilesional_distance+1, :, :);
        end
        
        FCD_IIa = strncmp(histo_pat, 'IIA', 6)';
        FCD_IIb = strncmp(histo_pat, 'IIB', 6)';
        
        p_value_set = [];
        
        % fMRI ALFF
        for fMRI_ALFF = 1
            
            feature_idx = 1; surface_idx = [1];
            x_lim = [0 perilesional_distance+2];
            y_lim = [-3 3];
            sig = -2;
            
            exclude_IIA = []; IDX_IIA = find(strcmp(histo_pat, 'IIA')); IDX_IIA(exclude_IIA) = [];
            exclude_IIB = []; IDX_IIB = find(strcmp(histo_pat, 'IIB')); IDX_IIB(exclude_IIB) = [];
            feature_type_IIA = mean(pl_feature_set(feature_idx).data(:, IDX_IIA, surface_idx), 3);
            feature_type_IIB = mean(pl_feature_set(feature_idx).data(:, IDX_IIB, surface_idx), 3);
            feature_type_IIA_ALFF = feature_type_IIA;
            feature_type_IIB_ALFF = feature_type_IIB;
            
            p_type_IIa = [];
            p_type_IIb = [];
            p_type_IIa_vs_IIb = [];
            for i = 1 : perilesional_distance+1
                mean_z_lesion_temp_a = mean(feature_type_IIA(i, :), 2);
                std_z_lesion_temp_a  = std(feature_type_IIA(i, :));
                pooled_std = sqrt(power(std_z_lesion_temp_a, 2)/sum(FCD_IIa)+1/size(case_num_cont, 1));
                t = (mean_z_lesion_temp_a - 0)/pooled_std;
                df = power(power(std_z_lesion_temp_a, 2)/sum(FCD_IIa)+1/size(case_num_cont, 1), 2)/(power(power(std_z_lesion_temp_a, 2)/sum(FCD_IIa),2)/(sum(FCD_IIa)-1) + power(1/size(case_num_cont, 1), 2)/(size(case_num_cont, 1)-1));
                p = tcdf(t, df);
                p_type_IIa = [ p_type_IIa min(p, 1-p) ];
                
                mean_z_lesion_temp_b = mean(feature_type_IIB(i, :), 2);
                std_z_lesion_temp_b  = std(feature_type_IIB(i, :));
                pooled_std = sqrt(power(std_z_lesion_temp_b, 2)/sum(FCD_IIb)+1/size(case_num_cont, 1));
                t = (mean_z_lesion_temp_b - 0)/pooled_std;
                df = power(power(std_z_lesion_temp_b, 2)/sum(FCD_IIb)+1/size(case_num_cont, 1), 2)/(power(power(std_z_lesion_temp_b, 2)/sum(FCD_IIb),2)/(sum(FCD_IIb)-1) + power(1/size(case_num_cont, 1), 2)/(size(case_num_cont, 1)-1));
                p = tcdf(t, df);
                p_type_IIb = [ p_type_IIb min(p, 1-p) ];
                
                pooled_std = sqrt(power(std_z_lesion_temp_a, 2)/sum(FCD_IIa)+power(std_z_lesion_temp_b, 2)/sum(FCD_IIb));
                t = (mean_z_lesion_temp_a - mean_z_lesion_temp_b)/pooled_std;
                df = power(power(std_z_lesion_temp_a, 2)/sum(FCD_IIa)+power(std_z_lesion_temp_b, 2)/sum(FCD_IIb), 2)/...
                    (power(power(std_z_lesion_temp_a, 2)/sum(FCD_IIa),2)/(sum(FCD_IIa)-1) + power(power(std_z_lesion_temp_b, 2)/sum(FCD_IIb),2)/(sum(FCD_IIb)-1));
                p = tcdf(t, df);
                p_type_IIa_vs_IIb = [ p_type_IIa_vs_IIb min(p, 1-p) ];
            end
            p_value_set = [ p_value_set p_type_IIa p_type_IIb p_type_IIa_vs_IIb];
            p_value_set_ALFF = [ p_type_IIa p_type_IIb p_type_IIa_vs_IIb];
            
            interval = 4;
            
            if(visualization)
                visualization_graph_perilesional_profile(perilesional_distance, p_type_IIa, p_type_IIb, feature_type_IIA, feature_type_IIB, x_lim, y_lim, sig, FDR_threshold_final_perilesion_profiling, interval)
                if(save_file)
                    export_fig([OUTPATH '/38_perilesional_profiling_ALFF_sm_' num2str(Kernel) ], '-m3', '-png' ); close(gcf);
                end
            end
            
            %% correlation to lesion size
            sig_range = [2:12];

            [r p] = corr(sum(lesion_label_data, 2), mean(mean(pl_feature_set(feature_idx).data(sig_range, :, surface_idx), 3), 1)')                        
            [r p] = corr(log(sum(lesion_label_data(IDX_IIA,:), 2)), mean(mean(pl_feature_set(feature_idx).data(sig_range, IDX_IIA, surface_idx), 3), 1)')
            [r p] = corr(sum(lesion_label_data(IDX_IIB,:), 2), mean(mean(pl_feature_set(feature_idx).data(sig_range, IDX_IIB, surface_idx), 3), 1)')
             
        end
        
        % fMRI ReHo
        for fMRI_ReHo = 1
            
            feature_idx = 2; surface_idx = [1];
            x_lim = [0 perilesional_distance+2];
            y_lim = [-3 3];
            sig = -2;
            
            exclude_IIA = []; IDX_IIA = find(strcmp(histo_pat, 'IIA')); IDX_IIA(exclude_IIA) = [];
            exclude_IIB = []; IDX_IIB = find(strcmp(histo_pat, 'IIB')); IDX_IIB(exclude_IIB) = [];
            feature_type_IIA = mean(pl_feature_set(feature_idx).data(:, IDX_IIA, surface_idx), 3);
            feature_type_IIB = mean(pl_feature_set(feature_idx).data(:, IDX_IIB, surface_idx), 3);
            feature_type_IIA_ReHo = feature_type_IIA;
            feature_type_IIB_ReHo = feature_type_IIB;
            
            p_type_IIa = [];
            p_type_IIb = [];
            p_type_IIa_vs_IIb = [];
            for i = 1 : perilesional_distance+1
                mean_z_lesion_temp_a = mean(feature_type_IIA(i, :), 2);
                std_z_lesion_temp_a  = std(feature_type_IIA(i, :));
                pooled_std = sqrt(power(std_z_lesion_temp_a, 2)/sum(FCD_IIa)+1/size(case_num_cont, 1));
                t = (mean_z_lesion_temp_a - 0)/pooled_std;
                df = power(power(std_z_lesion_temp_a, 2)/sum(FCD_IIa)+1/size(case_num_cont, 1), 2)/(power(power(std_z_lesion_temp_a, 2)/sum(FCD_IIa),2)/(sum(FCD_IIa)-1) + power(1/size(case_num_cont, 1), 2)/(size(case_num_cont, 1)-1));
                p = tcdf(t, df);
                p_type_IIa = [ p_type_IIa min(p, 1-p) ];
                
                mean_z_lesion_temp_b = mean(feature_type_IIB(i, :), 2);
                std_z_lesion_temp_b  = std(feature_type_IIB(i, :));
                pooled_std = sqrt(power(std_z_lesion_temp_b, 2)/sum(FCD_IIb)+1/size(case_num_cont, 1));
                t = (mean_z_lesion_temp_b - 0)/pooled_std;
                df = power(power(std_z_lesion_temp_b, 2)/sum(FCD_IIb)+1/size(case_num_cont, 1), 2)/(power(power(std_z_lesion_temp_b, 2)/sum(FCD_IIb),2)/(sum(FCD_IIb)-1) + power(1/size(case_num_cont, 1), 2)/(size(case_num_cont, 1)-1));
                p = tcdf(t, df);
                p_type_IIb = [ p_type_IIb min(p, 1-p) ];
                
                pooled_std = sqrt(power(std_z_lesion_temp_a, 2)/sum(FCD_IIa)+power(std_z_lesion_temp_b, 2)/sum(FCD_IIb));
                t = (mean_z_lesion_temp_a - mean_z_lesion_temp_b)/pooled_std;
                df = power(power(std_z_lesion_temp_a, 2)/sum(FCD_IIa)+power(std_z_lesion_temp_b, 2)/sum(FCD_IIb), 2)/...
                    (power(power(std_z_lesion_temp_a, 2)/sum(FCD_IIa),2)/(sum(FCD_IIa)-1) + power(power(std_z_lesion_temp_b, 2)/sum(FCD_IIb),2)/(sum(FCD_IIb)-1));
                p = tcdf(t, df);
                p_type_IIa_vs_IIb = [ p_type_IIa_vs_IIb min(p, 1-p) ];
            end
            p_value_set = [ p_value_set p_type_IIa p_type_IIb p_type_IIa_vs_IIb];
            p_value_set_ReHo = [ p_type_IIa p_type_IIb p_type_IIa_vs_IIb];
            interval = 4;
            
            if(visualization)
                visualization_graph_perilesional_profile(perilesional_distance, p_type_IIa, p_type_IIb, feature_type_IIA, feature_type_IIB, x_lim, y_lim, sig, FDR_threshold_final_perilesion_profiling, interval)
                if(save_file)
                    export_fig([OUTPATH '/39_perilesional_profiling_ReHo_sm_' num2str(Kernel) ], '-m3', '-png'); close(gcf);
                end
            end
            
            %% correlation to lesion size
            sig_range = [2:12];

            [r p] = corr(sum(lesion_label_data, 2), mean(mean(pl_feature_set(feature_idx).data(sig_range, :, surface_idx), 3), 1)')                        
            [r p] = corr(log(sum(lesion_label_data(IDX_IIA,:), 2)), mean(mean(pl_feature_set(feature_idx).data(sig_range, IDX_IIA, surface_idx), 3), 1)')
            [r p] = corr(sum(lesion_label_data(IDX_IIB,:), 2), mean(mean(pl_feature_set(feature_idx).data(sig_range, IDX_IIB, surface_idx), 3), 1)')
            
        end
        
    end
    
    if control_comp == 1
        
        for i = 1 : size(pl_feature_set_org, 2)
            temp = [];
            if(feature_parameters{i, 2} > 1)
                for c = 1 : size(case_num_cont, 1)
                    temp(:, c, :) = mean(squeeze(pl_feature_set_org(i).data(1:perilesional_distance+1, c, :, :)), 3);
                end
            elseif(feature_parameters{i, 2} == 1)
                for c = 1 : size(case_num_cont, 1)
                    temp(:, c, :) = mean(squeeze(pl_feature_set_org(i).data(1:perilesional_distance+1, c, :, :)), 2);
                end
            end
            
            pl_feature_set(i).data = squeeze(temp);
        end        
        
        p_value_set = [];
        
        % fMRI ALFF
        for fMRI_ALFF = 1
                        
            feature_idx = 1; surface_idx = [1];
            x_lim = [0 perilesional_distance+2];
            y_lim = [-2 2];
            sig = -2;            
           
            feature_type = mean(pl_feature_set(feature_idx).data(:, :, surface_idx), 3);                        
            
            if(visualization)
                visualization_graph_perilesional_profile_control(perilesional_distance, feature_type, x_lim, y_lim);
            end
            
        end
        
        % fMRI ReHo
        for fMRI_ReHo = 1
            
            feature_idx = 2; surface_idx = [1];
            x_lim = [0 perilesional_distance+2];
            y_lim = [-2 2];
            sig = -2;            
           
            feature_type = mean(pl_feature_set(feature_idx).data(:, :, surface_idx), 3);                        
            
            if(visualization)
                visualization_graph_perilesional_profile_control(perilesional_distance, feature_type, x_lim, y_lim);
            end
            
        end
        
    end
    
end