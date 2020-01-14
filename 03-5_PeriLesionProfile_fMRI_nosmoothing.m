clear all
close all
addpath('/data/noel/noel2/local/matlab/surfstat_chicago');

%% Read pre-calculated zscore of T1 and FLAIR in the lesion
for read_matfiles = 1
    
    load('zscore_database_perilesional_fMRI_nosmoothing.mat');
    load('individual_surface_multimodal.mat'); % geodesic distance computed based on original surface
                                               % and modified using small sphere kernel to solve jammed WM surfaces
    load([ 'p_value_set_perilesion_profile_T1_FLAIR.mat' ]);
    load([ 'p_value_set_perilesion_profile_DTI.mat' ]);
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
    pl_feature_read_again = 1;  % 0 = read presaved files | 1 = compute features again

    
    fid = fopen(Cases_cont);
    demo = textscan(fid, '%s%d%s%s%s', 'Delimiter', ',', 'CollectOutput', 1);
    case_num_cont = demo{1};
    age_cont = demo{2};
    gender_cont = demo{3}(:, 1);
    
    fclose(fid);
    
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
    
    % remove 061 075 079 from geoDist_dti
    [C, ia, ib] = intersect(case_num_pat_temp, case_num_pat);

    geoDist_fMRI     = geoDist_fMRI(ia);
    geoDist_fMRI_new = geoDist_fMRI_new(ia);
    
    surf_template = SurfStatReadSurf({'surf_reg_model_left.obj', 'surf_reg_model_right.obj'});
    mask_template = SurfStatMaskCut(surf_template);
    
    FDR_threshold_final_perilesion_profiling = FDR([p_value_perilesion_T1_FLAIR p_value_perilesion_DTI], 0.05); %% 0.0025
    
end

%% Compute the peri-lesional features along the cortical mantle
for feature_extraction = 1
    
   
    %% mean_z_lesion
    %% 1   = ALFF of fMRI              -> 1   p1
    %% 2   = ReHo of fMRI              -> 2   p2    
    
    feature_table = [];
    feature_parameters = [
        { 'fMRI_ALFF_z', 1, 1 }; %% 01 ALFF 1
        { 'fMRI_ReHo_z', 1, 1 }; %% 02 ReHO 2
        ];
    
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
    
    perilesional_distance_org = perilesional_distance;
    pl_feature_set_org = pl_feature_set;

end

%% Start perilesional analyses
for perilesional_analyses = 1
    
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
            p_type_IIa = [ p_type_IIa 2*min(p, 1-p) ];
            
            mean_z_lesion_temp_b = mean(feature_type_IIB(i, :), 2);
            std_z_lesion_temp_b  = std(feature_type_IIB(i, :));
            pooled_std = sqrt(power(std_z_lesion_temp_b, 2)/sum(FCD_IIb)+1/size(case_num_cont, 1));
            t = (mean_z_lesion_temp_b - 0)/pooled_std;
            df = power(power(std_z_lesion_temp_b, 2)/sum(FCD_IIb)+1/size(case_num_cont, 1), 2)/(power(power(std_z_lesion_temp_b, 2)/sum(FCD_IIb),2)/(sum(FCD_IIb)-1) + power(1/size(case_num_cont, 1), 2)/(size(case_num_cont, 1)-1));
            p = tcdf(t, df);
            p_type_IIb = [ p_type_IIb 2*min(p, 1-p) ];
            
            pooled_std = sqrt(power(std_z_lesion_temp_a, 2)/sum(FCD_IIa)+power(std_z_lesion_temp_b, 2)/sum(FCD_IIb));
            t = (mean_z_lesion_temp_a - mean_z_lesion_temp_b)/pooled_std;
            df = power(power(std_z_lesion_temp_a, 2)/sum(FCD_IIa)+power(std_z_lesion_temp_b, 2)/sum(FCD_IIb), 2)/...
                (power(power(std_z_lesion_temp_a, 2)/sum(FCD_IIa),2)/(sum(FCD_IIa)-1) + power(power(std_z_lesion_temp_b, 2)/sum(FCD_IIb),2)/(sum(FCD_IIb)-1));
            p = tcdf(t, df);
            p_type_IIa_vs_IIb = [ p_type_IIa_vs_IIb 2*min(p, 1-p) ];
        end
        p_value_set = [ p_value_set p_type_IIa p_type_IIb p_type_IIa_vs_IIb];
        interval = 4;
        
        if(visualization)
            visualization_graph_perilesional_profile(perilesional_distance, p_type_IIa, p_type_IIb, feature_type_IIA, feature_type_IIB, x_lim, y_lim, sig, FDR_threshold_final_perilesion_profiling, interval)
            if(save_file)
                export_fig([OUTPATH '/40_perilesional_profiling_ALFF_nosmoothing' ], '-m3', '-png'); close(gcf);
            end
        end
        
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
            p_type_IIa = [ p_type_IIa 2*min(p, 1-p) ];
            
            mean_z_lesion_temp_b = mean(feature_type_IIB(i, :), 2);
            std_z_lesion_temp_b  = std(feature_type_IIB(i, :));
            pooled_std = sqrt(power(std_z_lesion_temp_b, 2)/sum(FCD_IIb)+1/size(case_num_cont, 1));
            t = (mean_z_lesion_temp_b - 0)/pooled_std;
            df = power(power(std_z_lesion_temp_b, 2)/sum(FCD_IIb)+1/size(case_num_cont, 1), 2)/(power(power(std_z_lesion_temp_b, 2)/sum(FCD_IIb),2)/(sum(FCD_IIb)-1) + power(1/size(case_num_cont, 1), 2)/(size(case_num_cont, 1)-1));
            p = tcdf(t, df);
            p_type_IIb = [ p_type_IIb 2*min(p, 1-p) ];
            
            pooled_std = sqrt(power(std_z_lesion_temp_a, 2)/sum(FCD_IIa)+power(std_z_lesion_temp_b, 2)/sum(FCD_IIb));
            t = (mean_z_lesion_temp_a - mean_z_lesion_temp_b)/pooled_std;
            df = power(power(std_z_lesion_temp_a, 2)/sum(FCD_IIa)+power(std_z_lesion_temp_b, 2)/sum(FCD_IIb), 2)/...
                (power(power(std_z_lesion_temp_a, 2)/sum(FCD_IIa),2)/(sum(FCD_IIa)-1) + power(power(std_z_lesion_temp_b, 2)/sum(FCD_IIb),2)/(sum(FCD_IIb)-1));
            p = tcdf(t, df);
            p_type_IIa_vs_IIb = [ p_type_IIa_vs_IIb 2*min(p, 1-p) ];
        end
        p_value_set = [ p_value_set p_type_IIa p_type_IIb p_type_IIa_vs_IIb];
        interval = 4;
        
        if(visualization)
            visualization_graph_perilesional_profile(perilesional_distance, p_type_IIa, p_type_IIb, feature_type_IIA, feature_type_IIB, x_lim, y_lim, sig, FDR_threshold_final_perilesion_profiling, interval)
            if(save_file)
                export_fig([OUTPATH '/41_perilesional_profiling_ReHo_nosmoothing' ], '-m3', '-png'); close(gcf);
            end
        end
        
    end
    
end