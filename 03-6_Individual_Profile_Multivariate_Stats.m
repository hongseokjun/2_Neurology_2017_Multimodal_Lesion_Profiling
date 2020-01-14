clear all
close all
addpath('/data/noel/noel2/local/matlab/surfstat_chicago');

%% Read pre-calculated zscore of T1 and FLAIR in the lesion
for read_matfiles = 1
    
    load('zscore_database_T1_FLAIR.mat');
    load('zscore_database_T1_FLAIR_cont.mat');
    mean_z_lesion_T1_FLAIR      = mean_z_lesion;
    mean_z_lesion_T1_FLAIR_cont = mean_z_lesion_cont;
    
    load('zscore_database_DTI.mat');
    load('zscore_database_DTI_cont.mat');  
    mean_z_lesion_DTI      = mean_z_lesion;
    mean_z_lesion_DTI_cont = mean_z_lesion_cont;
    
    load('zscore_database_fMRI.mat');
    load('zscore_database_fMRI_cont.mat');  
    mean_z_lesion_fMRI      = mean_z_lesion;
    mean_z_lesion_fMRI_cont = mean_z_lesion_cont;
    
    load('/local_raid/seokjun/01_project/04_IntracorticalAnalysis/01_analysis/pipeline/postprocessing_zscore_visualization/colormap_noel1.mat');
    
end

%% Read demograpic data and set up some parameters and variables
for read_demodata_setup_params = 1
    
    OUTPATH = '/local_raid/seokjun/01_project/04_IntracorticalAnalysis/03_result/prediction_histology_by_MRI/';
    Cases_cont                = '/local_raid/seokjun/01_project/04_IntracorticalAnalysis/01_analysis/Demographic_data_control.txt';
    Cases_cont_fMRI           = '/local_raid/seokjun/01_project/04_IntracorticalAnalysis/01_analysis/Demographic_data_control_fMRI.txt';
    Group_cont                = 'control';
    Prefix_cont               = 'TLE';
    
    Cases_pat                 = '/local_raid/seokjun/01_project/04_IntracorticalAnalysis/01_analysis/Demographic_data_FCD.txt';
    Cases_pat_dti             = '/local_raid/seokjun/01_project/04_IntracorticalAnalysis/01_analysis/Demographic_data_FCD_DTI.txt';
    Cases_pat_fMRI            = '/local_raid/seokjun/01_project/04_IntracorticalAnalysis/01_analysis/Demographic_data_FCD_fMRI.txt';
    Group_pat                 = 'FCD';
    Prefix_pat                = 'mcd';
    
    Left_Right                = 'both';
    NumIntSurf                = 3;
    NumSubSurf                = 5;
    NumMesh                   = 81920;
    SamplingSpace             = 'native';
    img_contrast              = {'t1', 'flair', 'dti'};
    Kernel                    = 5;
    Parametric                = 'quadratic';
    average_surface_dir       = '/local_raid/seokjun/01_project/04_IntracorticalAnalysis/03_result/average_surfaces/';
    
    visualization = 1;
    save_fig = 0;
    load('/local_raid/seokjun/01_project/04_IntracorticalAnalysis/01_analysis/pipeline/postprocessing_zscore_visualization/colormap_noel1.mat');
    
    %% read controls data
    fid = fopen(Cases_cont);
    demo = textscan(fid, '%s%f%s', 'Delimiter', ',', 'CollectOutput', 1);
    case_num_cont = demo{1};
    age_cont = demo{2};
    gender_cont = demo{3};
    fclose(fid);
    
    %% read controls who have fMRI data
    fid = fopen(Cases_cont_fMRI);
    demo = textscan(fid, '%s%f%s%s%f', 'Delimiter', ',', 'CollectOutput', 1);
    case_num_cont_fMRI = demo{1}; 
    age_cont_fMRI = demo{2}; 
    gender_cont_fMRI = demo{3}(:, 1); 
    group_cont_fMRI = demo{3}(:, 2); 
    fclose(fid);
    
    excluded_controls = { '306_1', '313_1', '322_1',  '342_1', '343_1'};
    [C, ia, ib] = intersect(case_num_cont, excluded_controls);
    case_num_cont(ia) = [];
    age_cont(ia) = [];
    gender_cont(ia) = [];    
    
    excluded_controls = { '306_1', '313_1', '322_2',  '342_1', '343_1'};
    [C, ia, ib] = intersect(case_num_cont_fMRI, excluded_controls);
    case_num_cont_fMRI(ia) = [];
    age_cont_fMRI(ia) = [];
    gender_cont_fMRI(ia) = [];
    group_cont_fMRI(ia) = [];
    mean_z_lesion_fMRI_cont(:, :, ia) = [];
    
    %% read patients data (t1, flair)
    fid = fopen(Cases_pat);
    demo = textscan(fid, '%s%f%s%d%s%d%d%s%s', 'Delimiter', ',', 'CollectOutput', 1);
    case_num_pat = demo{1};
    age_pat = demo{2};
    gender_pat = demo{3};
    lesion_volume = demo{4};
    histo_type = demo{5};
    initial = demo{6}(:, 1);
    transmantle = demo{6}(:, 2);
    location = demo{7}(:, 1);
    seizure_lateralization = demo{7}(:, 2);
    fclose(fid);
    
    %% read patients data (dti)
    fid = fopen(Cases_pat_dti);
    demo = textscan(fid, '%s%f%s%d%s%d%d%s%s', 'Delimiter', ',', 'CollectOutput', 1);
    case_num_pat_dti = demo{1};
    age_pat_dti = demo{2};
    gender_pat_dti = demo{3};
    lesion_volume_dti = demo{4};
    histo_type_dti = demo{5};
    initial_dti = demo{6}(:, 1);
    transmantle_dti = demo{6}(:, 2);
    location_dti = demo{7}(:, 1);
    seizure_lateralization_dti = demo{7}(:, 2);
    fclose(fid);
        
    %% read patients data (fMRI)
    fid = fopen(Cases_pat_fMRI);
    demo = textscan(fid, '%s%d%s%s%s', 'Delimiter', ',', 'CollectOutput', 1);
    case_num_pat_fMRI = demo{1}; 
    age_pat_fMRI      = demo{2};
    gender_pat_fMRI   = demo{3}(:, 1);
    histo_pat_fMRI    = demo{3}(:, 3);
    histo_type_fMRI   = histo_pat_fMRI;
    histo_type_fMRI(strcmp(histo_pat_fMRI, 'IIB')) = {'FCDIIb'};
    histo_type_fMRI(strcmp(histo_pat_fMRI, 'IIA')) = {'FCDIIa'};
    
    %% select only patients who have dti data        
    case_num_pat{strcmp(case_num_pat, '080_1')} = '080';
    case_num_pat_dti{strcmp(case_num_pat_dti, '080_2')} = '080';
    case_num_pat_fMRI{strcmp(case_num_pat_fMRI, '080_2')} = '080';
        
end

%% Select features to be included in the analyses
for feature_selection = 1    
    
    [C, ia1, ib1] = intersect(case_num_pat, case_num_pat_dti);
    [C, ia2, ib2] = intersect(case_num_pat, case_num_pat_fMRI);
    [C, ia3, ib3] = intersect(case_num_pat_dti, case_num_pat_fMRI);
    [C, ia4, ib4] = intersect(case_num_cont, case_num_cont_fMRI);    
    
    FCD_IIa_org = strncmp(histo_type, 'FCDIIa', 6)';
    FCD_IIb_org = strncmp(histo_type, 'FCDIIb', 6)';
    FCD_IIa_dti = strncmp(histo_type_dti, 'FCDIIa', 6)';
    FCD_IIb_dti = strncmp(histo_type_dti, 'FCDIIb', 6)';
    FCD_IIa_fMRI = strncmp(histo_type_fMRI, 'FCDIIa', 6)';
    FCD_IIb_fMRI = strncmp(histo_type_fMRI, 'FCDIIb', 6)';
    
    mean_z_lesion_           = [ mean_z_lesion_T1_FLAIR ];
    mean_z_lesion_DTI_       = [ mean_z_lesion_T1_FLAIR(:, ia1); mean_z_lesion_DTI ];
    mean_z_lesion_fMRI_      = [ mean_z_lesion_T1_FLAIR(:, ia2); mean_z_lesion_DTI(:, ia3); mean_z_lesion_fMRI ];
    mean_z_lesion_(mean_z_lesion_>100) = 10;
    mean_z_lesion_DTI_(mean_z_lesion_DTI_>100) = 10;
    mean_z_lesion_fMRI_(mean_z_lesion_fMRI_>100) = 10;
    
    mean_z_lesion_cont_           = cat(1, mean_z_lesion_T1_FLAIR_cont);
    mean_z_lesion_DTI_cont_       = cat(1, mean_z_lesion_T1_FLAIR_cont(:, ia1, :), mean_z_lesion_DTI_cont);
    mean_z_lesion_fMRI_cont_      = cat(1, mean_z_lesion_T1_FLAIR_cont(:, ia2, ia4), mean_z_lesion_DTI_cont(:, ia3, ia4), mean_z_lesion_fMRI_cont);
    
    %% mean_z_lesion
    %% 1:5   = RI of T1                -> 2:5      p1
    %% 6:10  = PG of T1                -> 7:9 16   p2
    %% 11:15 = TG of T1                -> 12:15    p3
    %% 16    = PG_gw of T1             -> 7:9 16   p2
    %% 17    = CT of T1                -> 17       p7
    %% 18    = MC of T1                -> 18       p8
    %% 19    = SD of T1                -> 19       p9
    %% 20:22 = RI of T1 subcortical    -> 20:22    p4
    %% 23:25 = PG of T1 subcortical    -> 23:25    p5
    %% 26:28 = TG of T1 subcortical    -> 26:28    p6
    %% 29:33 = RI of FLAIR             -> 30:33    p10
    %% 34:38 = PG of FLAIR             -> 35:37 44 p11
    %% 39:43 = TG of FLAIR             -> 40:43    p12
    %% 44    = PG_gw of FLAIR          -> 35:37 44 p11
    %% 45:47 = RI of FLAIR subcortical -> 45:47    p13
    %% 48:50 = PG of FLAIR subcortical -> 48:50    p14
    %% 51:53 = TG of FLAIR subcortical -> 51:53    p15    
    feature_set_idx_T1_FLAIR = { 2:5, [7:9 16], 12:15, 20:22, 17, 18, 19, 30:33, [35:37 44], 40:43, 45:47 };
    
    %% mean_z_lesion
    %% 1:5   = FA of DTI Intra         -> 3 5       p1
    %% 6:12  = FA of DTI Subcortical   -> 7 9 11    p2
    %% 13:17 = MD of DTI Intra         -> 15 17     p3
    %% 18:24 = MD of DTI Subcortical   -> 19 21 23  p4
    feature_set_idx_DTI      = { 3+size(mean_z_lesion_T1_FLAIR, 1), 5+size(mean_z_lesion_T1_FLAIR, 1), 7+size(mean_z_lesion_T1_FLAIR, 1), ...
                                 9+size(mean_z_lesion_T1_FLAIR, 1), 11+size(mean_z_lesion_T1_FLAIR, 1), 15+size(mean_z_lesion_T1_FLAIR, 1), ...
                                 17+size(mean_z_lesion_T1_FLAIR, 1), 19+size(mean_z_lesion_T1_FLAIR, 1), 21+size(mean_z_lesion_T1_FLAIR, 1), 23+size(mean_z_lesion_T1_FLAIR, 1) };
    
    %% mean_z_lesion
    %% 1     = ALFF                    -> 1       p1
    %% 2     = ReHo                    -> 2       p2
    %% 3     = DC                      -> 3       p3
    %% 4     = BC                      -> 4       p4
    feature_set_idx_fMRI     = { 1+size(mean_z_lesion_T1_FLAIR, 1)+size(mean_z_lesion_DTI, 1), ...
                                 2+size(mean_z_lesion_T1_FLAIR, 1)+size(mean_z_lesion_DTI, 1), ...
                                 3+size(mean_z_lesion_T1_FLAIR, 1)+size(mean_z_lesion_DTI, 1), ...
                                 4+size(mean_z_lesion_T1_FLAIR, 1)+size(mean_z_lesion_DTI, 1)};
end

%% Start individual analysis
for individual_analyses = 1
    
    %% mean_z_lesion
    %% 1:5   = RI of T1                -> 2:5      p1
    %% 6:10  = PG of T1                -> 7:9 16   p2
    %% 11:15 = TG of T1                -> 12:15    p3
    %% 16    = PG_gw of T1             -> 7:9 16   
    %% 17    = CT of T1                -> 17       p5
    %% 18    = MC of T1                -> 18       p6
    %% 19    = SD of T1                -> 19       p7
    %% 20:22 = RI of T1 subcortical    -> 20:22    p4
    %% 23:25 = PG of T1 subcortical    -> 23:25    
    %% 26:28 = TG of T1 subcortical    -> 26:28    
    %% 29:33 = RI of FLAIR             -> 30:33    p8
    %% 34:38 = PG of FLAIR             -> 35:37 44 p9
    %% 39:43 = TG of FLAIR             -> 40:43    p10
    %% 44    = PG_gw of FLAIR          -> 35:37 44 
    %% 45:47 = RI of FLAIR subcortical -> 45:47    p11
    %% 48:50 = PG of FLAIR subcortical -> 48:50    
    %% 51:53 = TG of FLAIR subcortical -> 51:53            
    %% 54:58 = FA of DTI Intra         -> 56 58    p12
    %% 59:65 = FA of DTI Subcortical   -> 60 62 64 p14
    %% 66:70 = MD of DTI Intra         -> 68 70    p13
    %% 71:77 = MD of DTI Subcortical   -> 72 74 76 p15
      
    FCD_IIa_T1_FLAIR    = strncmp(histo_type,           'FCDIIa', 6)';
    FCD_IIb_T1_FLAIR    = strncmp(histo_type,           'FCDIIb', 6)';
    FCDIIa_T1_FLAIR_idx = find(strcmp(histo_type,       'FCDIIa'));
    FCDIIb_T1_FLAIR_idx = find(strcmp(histo_type,       'FCDIIb'));
    
    FCD_IIa_DTI         = strncmp(histo_type_dti,       'FCDIIa', 6)';
    FCD_IIb_DTI         = strncmp(histo_type_dti,       'FCDIIb', 6)';
    FCDIIa_DTI_idx      = find(strcmp(histo_type_dti,   'FCDIIa'));
    FCDIIb_DTI_idx      = find(strcmp(histo_type_dti,   'FCDIIb'));
    
    FCD_IIa_fMRI        = strncmp(histo_type_fMRI,      'FCDIIa', 6)';
    FCD_IIb_fMRI        = strncmp(histo_type_fMRI,      'FCDIIb', 6)';
    FCDIIa_fMRI_idx     = find(strcmp(histo_type_fMRI,  'FCDIIa'));
    FCDIIb_fMRI_idx     = find(strcmp(histo_type_fMRI,  'FCDIIb'));
    
    for by_modalities = 1
        
                          %   1       2         3        4      5   6   7    8         9           10      11        12      13       14       15     16  17
        feature_set     = { [2:5], [7:9 16], [12:15], [20:22], 17, 18, 19, [30:33], [35:37 44], [40:43], [45:47], [56 58], [68 70], [60 62], [72 74], 78, 79 };
        feature_set_str = { 'T1_cortRI', 'T1_PG', 'T1_TG', 'T1_subcRI', 'CT', 'Curv', 'SD', 'FL_cortRI', 'FL_PG', 'FL_TG', 'FL_subcRI', 'DTI_cortFA', 'DTI_cortMD', 'DTI_subcFA', 'DTI_subcMD' };
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % IIa_individual_analysis
        for IIa_individual_analysis = 1
            
            %% T1/FLAIR
            for T1_FLAIR = 1
                threshold_sign = [ 1 -1 -1 -1 1 1 1 1 -1 -1 1 -1 1 -1 1 -1 1];
                IIa_threshold = 1.96;
                
                % each feature
                feature_FCD_IIa_idx = [ 9 ];
                
                feature_FCD_IIa_set = feature_set(feature_FCD_IIa_idx);
                IIa_threshold_sign = threshold_sign(feature_FCD_IIa_idx);
                IIa_threshold_set = IIa_threshold_sign*IIa_threshold;
                feature_FCD_IIa         = [];
                feature_FCD_IIa_cont    = [];
                
                for i = 1 : length(feature_FCD_IIa_idx)
                    
                    feature_FCD_IIa             = [ feature_FCD_IIa         mean(mean_z_lesion_(feature_FCD_IIa_set{i}, FCD_IIa_T1_FLAIR), 1)' ];
                    if(length(feature_FCD_IIa_set{i}) > 1)
                        feature_FCD_IIa_cont    = [ feature_FCD_IIa_cont    mean(squeeze(mean(mean_z_lesion_cont_(feature_FCD_IIa_set{i}, FCD_IIa_T1_FLAIR, :), 2)), 1)' ];
                    else
                        feature_FCD_IIa_cont    = [ feature_FCD_IIa_cont    mean(squeeze(mean_z_lesion_cont_(feature_FCD_IIa_set{i}, FCD_IIa_T1_FLAIR, :)), 1)' ];
                    end
                    
                end
                
                abnorm = sum(sum(abs(feature_FCD_IIa)>IIa_threshold, 2)>=1);
                disp(['n=' num2str(abnorm), ' (' num2str(abnorm/sum(FCD_IIa_T1_FLAIR)) '%)']);
                
                %% multivariate features
                mahal_IIa      = mahal(feature_FCD_IIa, feature_FCD_IIa_cont);
                threshold_IIa  = mahal(IIa_threshold_set, feature_FCD_IIa_cont);
                disp(['n=' num2str(sum(mahal_IIa>=threshold_IIa)), ' (' num2str(sum(mahal_IIa>=threshold_IIa)/sum(FCD_IIa_T1_FLAIR)) '%)']);
                
                [a b] = max(mahal_IIa);
                case_num_pat_fMRI(FCDIIa_idx(b)) % highest md score: 065 071
                
                mahal_cont_IIa = [];
                for i = 1 : length(feature_FCD_IIa_cont)
                    feature_FCD_IIa_cont_temp = feature_FCD_IIa_cont;
                    feature_FCD_IIa_cont_temp(i, :) = [];
                    mahal_cont_IIa(i, :) = mahal(feature_FCD_IIa_cont(i, :), feature_FCD_IIa_cont_temp);
                end
                sum(mahal_cont_IIa>=threshold_IIa)/length(feature_FCD_IIa_cont)
            end
            
            %% DTI
            for DTI = 1
                threshold_sign = [ 1 -1 -1 -1 1 1 1 1 -1 -1 1 -1 1 -1 1 -1 1];
                IIa_threshold = 1.96;
                
                % each feature
                feature_FCD_IIa_idx = [ 14 15 ];
                
                feature_FCD_IIa_set = feature_set(feature_FCD_IIa_idx);
                IIa_threshold_sign = threshold_sign(feature_FCD_IIa_idx);
                IIa_threshold_set = IIa_threshold_sign*IIa_threshold;
                feature_FCD_IIa         = [];
                feature_FCD_IIa_cont    = [];
                
                for i = 1 : length(feature_FCD_IIa_idx)
                    
                    feature_FCD_IIa             = [ feature_FCD_IIa         mean(mean_z_lesion_DTI_(feature_FCD_IIa_set{i}, FCD_IIa_DTI), 1)' ];
                    if(length(feature_FCD_IIa_set{i}) > 1)
                        feature_FCD_IIa_cont    = [ feature_FCD_IIa_cont    mean(squeeze(mean(mean_z_lesion_DTI_cont_(feature_FCD_IIa_set{i}, FCD_IIa_DTI, :), 2)), 1)' ];
                    else
                        feature_FCD_IIa_cont    = [ feature_FCD_IIa_cont    mean(squeeze(mean_z_lesion_DTI_cont_(feature_FCD_IIa_set{i}, FCD_IIa_DTI, :)), 1)' ];
                    end
                    
                end
                
                abnorm = sum(sum(abs(feature_FCD_IIa)>IIa_threshold, 2)>=1);
                disp(['n=' num2str(abnorm), ' (' num2str(abnorm/sum(FCD_IIa_dti)) '%)']);
                
                %% multivariate features
                mahal_IIa      = mahal(feature_FCD_IIa, feature_FCD_IIa_cont);
                threshold_IIa  = mahal(IIa_threshold_set, feature_FCD_IIa_cont);
                disp(['n=' num2str(sum(mahal_IIa>=threshold_IIa)), ' (' num2str(sum(mahal_IIa>=threshold_IIa)/sum(FCD_IIa_DTI)) '%)']);
                
                [a b] = max(mahal_IIa);
                case_num_pat_fMRI(FCDIIa_idx(b)) % highest md score: 065 071
                
                mahal_cont_IIa = [];
                for i = 1 : length(feature_FCD_IIa_cont)
                    feature_FCD_IIa_cont_temp = feature_FCD_IIa_cont;
                    feature_FCD_IIa_cont_temp(i, :) = [];
                    mahal_cont_IIa(i, :) = mahal(feature_FCD_IIa_cont(i, :), feature_FCD_IIa_cont_temp);
                end
                sum(mahal_cont_IIa>=threshold_IIa)/length(feature_FCD_IIa_cont)
            end
            
            %% fMRI
            for fMRI = 1
                threshold_sign = [ 1 -1 -1 -1 1 1 1 1 -1 -1 1 -1 1 -1 1 -1 1];
                IIa_threshold = 1.96;
                
                % each feature
                feature_FCD_IIa_idx = [ 17 ];
                
                feature_FCD_IIa_set = feature_set(feature_FCD_IIa_idx);
                IIa_threshold_sign = threshold_sign(feature_FCD_IIa_idx);
                IIa_threshold_set = IIa_threshold_sign*IIa_threshold;
                feature_FCD_IIa         = [];
                feature_FCD_IIa_cont    = [];
                
                for i = 1 : length(feature_FCD_IIa_idx)
                    
                    feature_FCD_IIa             = [ feature_FCD_IIa         mean(mean_z_lesion_fMRI_(feature_FCD_IIa_set{i}, FCD_IIa_fMRI), 1)' ];
                    if(length(feature_FCD_IIa_set{i}) > 1)
                        feature_FCD_IIa_cont    = [ feature_FCD_IIa_cont    mean(squeeze(mean(mean_z_lesion_fMRI_cont_(feature_FCD_IIa_set{i}, FCD_IIa_fMRI, :), 2)), 1)' ];
                    else
                        feature_FCD_IIa_cont    = [ feature_FCD_IIa_cont    mean(squeeze(mean_z_lesion_fMRI_cont_(feature_FCD_IIa_set{i}, FCD_IIa_fMRI, :)), 1)' ];
                    end
                    
                end
                
                abnorm = sum(sum(abs(feature_FCD_IIa)>IIa_threshold, 2)>=1);
                disp(['n=' num2str(abnorm), ' (' num2str(abnorm/sum(FCD_IIa_fMRI)) '%)']);
                
                %% multivariate features
                mahal_IIa      = mahal(feature_FCD_IIa, feature_FCD_IIa_cont);
                threshold_IIa  = mahal(IIa_threshold_set, feature_FCD_IIa_cont);
                disp(['n=' num2str(sum(mahal_IIa>=threshold_IIa)), ' (' num2str(sum(mahal_IIa>=threshold_IIa)/sum(FCD_IIa_fMRI)) '%)']);
                
                [a b] = max(mahal_IIa);
                case_num_pat_fMRI(FCDIIa_idx(b)) % highest md score: 065 071
                
                mahal_cont_IIa = [];
                for i = 1 : length(feature_FCD_IIa_cont)
                    feature_FCD_IIa_cont_temp = feature_FCD_IIa_cont;
                    feature_FCD_IIa_cont_temp(i, :) = [];
                    mahal_cont_IIa(i, :) = mahal(feature_FCD_IIa_cont(i, :), feature_FCD_IIa_cont_temp);
                end
                sum(mahal_cont_IIa>=threshold_IIa)/length(feature_FCD_IIa_cont)
            end
            
            %% All
            for All = 1
                threshold_sign = [ 1 -1 -1 -1 1 1 1 1 -1 -1 1 -1 1 -1 1 -1 1];
                IIa_threshold = 1.5;
                
                feature_FCD_IIa_idx = [2 3 5 9 10 11 14];
                
                feature_FCD_IIa_set = feature_set(feature_FCD_IIa_idx);
                IIa_threshold_sign = threshold_sign(feature_FCD_IIa_idx);
                IIa_threshold_set = IIa_threshold_sign*IIa_threshold;
                feature_FCD_IIa         = [];
                feature_FCD_IIa_cont    = [];
                
                for i = 1 : length(feature_FCD_IIa_idx)
                    
                    feature_FCD_IIa             = [ feature_FCD_IIa         mean(mean_z_lesion_fMRI_(feature_FCD_IIa_set{i}, FCD_IIa_fMRI), 1)' ];
                    if(length(feature_FCD_IIa_set{i}) > 1)
                        feature_FCD_IIa_cont    = [ feature_FCD_IIa_cont    mean(squeeze(mean(mean_z_lesion_fMRI_cont_(feature_FCD_IIa_set{i}, FCD_IIa_fMRI, :), 2)), 1)' ];
                    else
                        feature_FCD_IIa_cont    = [ feature_FCD_IIa_cont    mean(squeeze(mean_z_lesion_fMRI_cont_(feature_FCD_IIa_set{i}, FCD_IIa_fMRI, :)), 1)' ];
                    end
                    
                end
                
                mahal_IIa      = mahal(feature_FCD_IIa, feature_FCD_IIa_cont);
                threshold_IIa  = mahal(IIa_threshold_set, feature_FCD_IIa_cont);
                disp(['n=' num2str(sum(mahal_IIa>=threshold_IIa)), ' (' num2str(sum(mahal_IIa>=threshold_IIa)/sum(FCD_IIa_fMRI)) '%)']);
                
                [a b] = max(mahal_IIa);
                case_num_pat_fMRI(FCDIIa_fMRI_idx(b)) % highest md score: 065 071
                
                mahal_cont_IIa = [];
                for i = 1 : length(feature_FCD_IIa_cont)
                    feature_FCD_IIa_cont_temp = feature_FCD_IIa_cont;
                    feature_FCD_IIa_cont_temp(i, :) = [];
                    mahal_cont_IIa(i, :) = mahal(feature_FCD_IIa_cont(i, :), feature_FCD_IIa_cont_temp);
                end
                sum(mahal_cont_IIa>=threshold_IIa)/length(feature_FCD_IIa_cont)
                
            end
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % IIb_individual_analysis
        for IIb_individual_analysis = 1
            
            %% T1/FLAIR
            for T1_FLAIR = 1
                
                threshold_sign = [ 1 -1 -1 -1 1 1 1 1 -1 -1 1 -1 1 -1 1 -1 1];
                IIb_threshold = 1.96;
                
                % each feature
                feature_FCD_IIb_idx = [ 7 ];
                
                feature_FCD_IIb_set = feature_set(feature_FCD_IIb_idx);
                IIb_threshold_sign = threshold_sign(feature_FCD_IIb_idx);
                IIb_threshold_set = IIb_threshold_sign*IIb_threshold;
                feature_FCD_IIb         = [];
                feature_FCD_IIb_cont    = [];
                
                for i = 1 : length(feature_FCD_IIb_idx)
                    
                    feature_FCD_IIb             = [ feature_FCD_IIb         mean(mean_z_lesion_(feature_FCD_IIb_set{i}, FCD_IIb_T1_FLAIR), 1)' ];
                    if(length(feature_FCD_IIb_set{i}) > 1)
                        feature_FCD_IIb_cont    = [ feature_FCD_IIb_cont    mean(squeeze(mean(mean_z_lesion_cont_(feature_FCD_IIb_set{i}, FCD_IIb_T1_FLAIR, :), 2)), 1)' ];
                    else
                        feature_FCD_IIb_cont    = [ feature_FCD_IIb_cont    mean(squeeze(mean_z_lesion_cont_(feature_FCD_IIb_set{i}, FCD_IIb_T1_FLAIR, :)), 1)' ];
                    end
                    
                end
                
                abnorm = sum(sum(abs(feature_FCD_IIb)>IIb_threshold, 2)>=1);
                disp(['n=' num2str(abnorm), ' (' num2str(abnorm/sum(FCD_IIb_T1_FLAIR)) '%)']);
                
                %% multivariate features
                mahal_IIb      = mahal(feature_FCD_IIb, feature_FCD_IIb_cont);
                threshold_IIb  = mahal(IIb_threshold_set, feature_FCD_IIb_cont);
                disp(['n=' num2str(sum(mahal_IIb>=threshold_IIb)), ' (' num2str(sum(mahal_IIb>=threshold_IIb)/sum(FCD_IIb_T1_FLAIR)) '%)']);
                
                [a b] = max(mahal_IIb);
                case_num_pat_fMRI(FCDIIb_idx(b)) % highest md score: 065 071
                
                mahal_cont_IIb = [];
                for i = 1 : length(feature_FCD_IIb_cont)
                    feature_FCD_IIb_cont_temp = feature_FCD_IIb_cont;
                    feature_FCD_IIb_cont_temp(i, :) = [];
                    mahal_cont_IIb(i, :) = mahal(feature_FCD_IIb_cont(i, :), feature_FCD_IIb_cont_temp);
                end
                sum(mahal_cont_IIb>=threshold_IIb)/length(feature_FCD_IIb_cont)
                
            end
            
            %% DTI
            for DTI = 1
                
                threshold_sign = [ 1 -1 -1 -1 1 1 1 1 -1 -1 1 -1 1 -1 1 -1 1];
                IIb_threshold = 1.96;
                
                % each feature
                feature_FCD_IIb_idx = [ 14 15 ];
                
                feature_FCD_IIb_set = feature_set(feature_FCD_IIb_idx);
                IIb_threshold_sign = threshold_sign(feature_FCD_IIb_idx);
                IIb_threshold_set = IIb_threshold_sign*IIb_threshold;
                feature_FCD_IIb         = [];
                feature_FCD_IIb_cont    = [];
                
                for i = 1 : length(feature_FCD_IIb_idx)
                    
                    feature_FCD_IIb             = [ feature_FCD_IIb         mean(mean_z_lesion_DTI_(feature_FCD_IIb_set{i}, FCD_IIb_DTI), 1)' ];
                    if(length(feature_FCD_IIb_set{i}) > 1)
                        feature_FCD_IIb_cont    = [ feature_FCD_IIb_cont    mean(squeeze(mean(mean_z_lesion_DTI_cont_(feature_FCD_IIb_set{i}, FCD_IIb_DTI, :), 2)), 1)' ];
                    else
                        feature_FCD_IIb_cont    = [ feature_FCD_IIb_cont    mean(squeeze(mean_z_lesion_DTI_cont_(feature_FCD_IIb_set{i}, FCD_IIb_DTI, :)), 1)' ];
                    end
                    
                end
                
                abnorm = sum(sum(abs(feature_FCD_IIb)>IIb_threshold, 2)>=1);
                disp(['n=' num2str(abnorm), ' (' num2str(abnorm/sum(FCD_IIb_dti)) '%)']);
                
                %% multivariate features
                mahal_IIb      = mahal(feature_FCD_IIb, feature_FCD_IIb_cont);
                threshold_IIb  = mahal(IIb_threshold_set, feature_FCD_IIb_cont);
                disp(['n=' num2str(sum(mahal_IIb>=threshold_IIb)), ' (' num2str(sum(mahal_IIb>=threshold_IIb)/sum(FCD_IIb_DTI)) '%)']);
                
                [a b] = max(mahal_IIb);
                case_num_pat_fMRI(FCDIIb_idx(b)) % highest md score: 065 071
                
                mahal_cont_IIb = [];
                for i = 1 : length(feature_FCD_IIb_cont)
                    feature_FCD_IIb_cont_temp = feature_FCD_IIb_cont;
                    feature_FCD_IIb_cont_temp(i, :) = [];
                    mahal_cont_IIb(i, :) = mahal(feature_FCD_IIb_cont(i, :), feature_FCD_IIb_cont_temp);
                end
                sum(mahal_cont_IIb>=threshold_IIb)/length(feature_FCD_IIb_cont)
                
            end
            
            %% fMRI
            for fMRI = 1
                
                threshold_sign = [ 1 -1 -1 -1 1 1 1 1 -1 -1 1 -1 1 -1 1 -1 1];
                IIb_threshold = 1.96;
                
                % each feature
                feature_FCD_IIb_idx = [ 17 ];
                
                feature_FCD_IIb_set = feature_set(feature_FCD_IIb_idx);
                IIb_threshold_sign = threshold_sign(feature_FCD_IIb_idx);
                IIb_threshold_set = IIb_threshold_sign*IIb_threshold;
                feature_FCD_IIb         = [];
                feature_FCD_IIb_cont    = [];
                
                for i = 1 : length(feature_FCD_IIb_idx)
                    
                    feature_FCD_IIb             = [ feature_FCD_IIb         mean(mean_z_lesion_fMRI_(feature_FCD_IIb_set{i}, FCD_IIb_fMRI), 1)' ];
                    if(length(feature_FCD_IIb_set{i}) > 1)
                        feature_FCD_IIb_cont    = [ feature_FCD_IIb_cont    mean(squeeze(mean(mean_z_lesion_fMRI_cont_(feature_FCD_IIb_set{i}, FCD_IIb_fMRI, :), 2)), 1)' ];
                    else
                        feature_FCD_IIb_cont    = [ feature_FCD_IIb_cont    mean(squeeze(mean_z_lesion_fMRI_cont_(feature_FCD_IIb_set{i}, FCD_IIb_fMRI, :)), 1)' ];
                    end
                    
                end
                
                abnorm = sum(sum(abs(feature_FCD_IIb)>IIb_threshold, 2)>=1);
                disp(['n=' num2str(abnorm), ' (' num2str(abnorm/sum(FCD_IIb_fMRI)) '%)']);
                
                %% multivariate features
                mahal_IIb      = mahal(feature_FCD_IIb, feature_FCD_IIb_cont);
                threshold_IIb  = mahal(IIb_threshold_set, feature_FCD_IIb_cont);
                disp(['n=' num2str(sum(mahal_IIb>=threshold_IIb)), ' (' num2str(sum(mahal_IIb>=threshold_IIb)/sum(FCD_IIb_fMRI)) '%)']);
                
                [a b] = max(mahal_IIb);
                case_num_pat_fMRI(FCDIIb_idx(b)) % highest md score: 065 071
                
                mahal_cont_IIb = [];
                for i = 1 : length(feature_FCD_IIb_cont)
                    feature_FCD_IIb_cont_temp = feature_FCD_IIb_cont;
                    feature_FCD_IIb_cont_temp(i, :) = [];
                    mahal_cont_IIb(i, :) = mahal(feature_FCD_IIb_cont(i, :), feature_FCD_IIb_cont_temp);
                end
                sum(mahal_cont_IIb>=threshold_IIb)/length(feature_FCD_IIb_cont)
                
            end
            
            %% All
            for All = 1
                threshold_sign = [ 1 -1 -1 -1 1 1 1 1 -1 -1 1 -1 1 -1 1 -1 1];
                IIb_threshold = 1.5;
                
                feature_FCD_IIb_idx = [2 3 4 5 7 8 9 11 15];
                
                feature_FCD_IIb_set = feature_set(feature_FCD_IIb_idx);
                IIb_threshold_sign = threshold_sign(feature_FCD_IIb_idx);
                IIb_threshold_set = IIb_threshold_sign*IIb_threshold;
                feature_FCD_IIb         = [];
                feature_FCD_IIb_cont    = [];
                
                for i = 1 : length(feature_FCD_IIb_idx)
                    
                    feature_FCD_IIb             = [ feature_FCD_IIb         mean(mean_z_lesion_fMRI_(feature_FCD_IIb_set{i}, FCD_IIb_fMRI), 1)' ];
                    if(length(feature_FCD_IIb_set{i}) > 1)
                        feature_FCD_IIb_cont    = [ feature_FCD_IIb_cont    mean(squeeze(mean(mean_z_lesion_fMRI_cont_(feature_FCD_IIb_set{i}, FCD_IIb_fMRI, :), 2)), 1)' ];
                    else
                        feature_FCD_IIb_cont    = [ feature_FCD_IIb_cont    mean(squeeze(mean_z_lesion_fMRI_cont_(feature_FCD_IIb_set{i}, FCD_IIb_fMRI, :)), 1)' ];
                    end
                    
                end
                
                mahal_IIb      = mahal(feature_FCD_IIb, feature_FCD_IIb_cont);
                threshold_IIb  = mahal(IIb_threshold_set, feature_FCD_IIb_cont);
                disp(['n=' num2str(sum(mahal_IIb>=threshold_IIb)), ' (' num2str(sum(mahal_IIb>=threshold_IIb)/sum(FCD_IIb_fMRI)) '%)']);
                
                [a b] = max(mahal_IIb);
                case_num_pat_fMRI(FCDIIb_fMRI_idx(b)) % highest md score: 065 071
                
                mahal_cont_IIb = [];
                for i = 1 : length(feature_FCD_IIb_cont)
                    feature_FCD_IIb_cont_temp = feature_FCD_IIb_cont;
                    feature_FCD_IIb_cont_temp(i, :) = [];
                    mahal_cont_IIb(i, :) = mahal(feature_FCD_IIb_cont(i, :), feature_FCD_IIb_cont_temp);
                end
                sum(mahal_cont_IIb>=threshold_IIb)/length(feature_FCD_IIb_cont)
                
            end
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % all_individual_analysis
        for all_individual_analysis = 1
            
            %% T1/FLAIR
            for T1_FLAIR = 1
                
                threshold_sign = [ 1 -1 -1 -1 1 1 1 1 -1 -1 1 -1 1 -1 1 -1 1];
                                
                % each feature                
                feature_idx = [ 1 2 3 4 8 9 10 11 ];
                threshold_set1 = 2.0;
                threshold_set2 = 0.0;
                feature_set_temp = feature_set(feature_idx);
                feature_set_str_temp = feature_set_str(feature_idx);
                threshold_sign = threshold_sign(feature_idx);

                threshold_set2_ = threshold_sign*threshold_set2;
                                
                                
                feature_FCD         = [];
                feature_FCD_cont    = [];
                
                for i = 1 : length(feature_set_temp)
                    
                    feature_FCD             = [ feature_FCD mean(mean_z_lesion_(feature_set_temp{i}, :), 1)' ];
                    if(length(feature_set_temp{i}) > 1)
                        feature_FCD_cont    = [ feature_FCD_cont    mean(squeeze(mean(mean_z_lesion_cont_(feature_set_temp{i}, :, :), 2)), 1)' ];
                    else
                        feature_FCD_cont    = [ feature_FCD_cont    mean(squeeze(mean_z_lesion_cont_(feature_set_temp{i}, :, :)), 1)' ];
                    end
                    
                end
                            
                %% multivariate features
                mahal_FCD      = mahal(feature_FCD, feature_FCD_cont); 
                
                threshold_FCD_set = [];
                x = 0 : 0.01 : 5;
                for i = x
                    
                    threshold_set_ = threshold_sign*i;
                    threshold_FCD_set(round(i/0.01+1))  = mahal(threshold_set_, feature_FCD_cont);
                    
                end
                p = polyfit(x, threshold_FCD_set ,2);
                y_sim = polyval(p, x);
                norm_mahal_FCD = sqrt(mahal_FCD/p(1));
                
                T = array2table([ norm_mahal_FCD feature_FCD ], 'VariableNames', [ 'mahal', feature_set_str_temp ], 'RowNames', case_num_pat)               
                
                
            end
            
            %% DTI
            for DTI = 1
                
                threshold_sign = [ 1 -1 -1 -1 1 1 1 1 -1 -1 1 -1 1 -1 1 -1 1];
                IIb_threshold = 1.96;
                
                % each feature
                feature_FCD_IIb_idx = [ 14 15 ];
                
                feature_FCD_IIb_set = feature_set(feature_FCD_IIb_idx);
                IIb_threshold_sign = threshold_sign(feature_FCD_IIb_idx);
                IIb_threshold_set = IIb_threshold_sign*IIb_threshold;
                feature_FCD_IIb         = [];
                feature_FCD_IIb_cont    = [];
                
                for i = 1 : length(feature_FCD_IIb_idx)
                    
                    feature_FCD_IIb             = [ feature_FCD_IIb         mean(mean_z_lesion_DTI_(feature_FCD_IIb_set{i}, FCD_IIb_DTI), 1)' ];
                    if(length(feature_FCD_IIb_set{i}) > 1)
                        feature_FCD_IIb_cont    = [ feature_FCD_IIb_cont    mean(squeeze(mean(mean_z_lesion_DTI_cont_(feature_FCD_IIb_set{i}, FCD_IIb_DTI, :), 2)), 1)' ];
                    else
                        feature_FCD_IIb_cont    = [ feature_FCD_IIb_cont    mean(squeeze(mean_z_lesion_DTI_cont_(feature_FCD_IIb_set{i}, FCD_IIb_DTI, :)), 1)' ];
                    end
                    
                end
                
                abnorm = sum(sum(abs(feature_FCD_IIb)>IIb_threshold, 2)>=1);
                disp(['n=' num2str(abnorm), ' (' num2str(abnorm/sum(FCD_IIb_dti)) '%)']);
                
                %% multivariate features
                mahal_IIb      = mahal(feature_FCD_IIb, feature_FCD_IIb_cont);
                threshold_IIb  = mahal(IIb_threshold_set, feature_FCD_IIb_cont);
                disp(['n=' num2str(sum(mahal_IIb>=threshold_IIb)), ' (' num2str(sum(mahal_IIb>=threshold_IIb)/sum(FCD_IIb_DTI)) '%)']);
                
                [a b] = max(mahal_IIb);
                case_num_pat_fMRI(FCDIIb_idx(b)) % highest md score: 065 071
                
                mahal_cont_IIb = [];
                for i = 1 : length(feature_FCD_IIb_cont)
                    feature_FCD_IIb_cont_temp = feature_FCD_IIb_cont;
                    feature_FCD_IIb_cont_temp(i, :) = [];
                    mahal_cont_IIb(i, :) = mahal(feature_FCD_IIb_cont(i, :), feature_FCD_IIb_cont_temp);
                end
                sum(mahal_cont_IIb>=threshold_IIb)/length(feature_FCD_IIb_cont)
                
            end
            
            %% fMRI
            for fMRI = 1
                
                threshold_sign = [ 1 -1 -1 -1 1 1 1 1 -1 -1 1 -1 1 -1 1 -1 1];
                IIb_threshold = 1.96;
                
                % each feature
                feature_FCD_IIb_idx = [ 17 ];
                
                feature_FCD_IIb_set = feature_set(feature_FCD_IIb_idx);
                IIb_threshold_sign = threshold_sign(feature_FCD_IIb_idx);
                IIb_threshold_set = IIb_threshold_sign*IIb_threshold;
                feature_FCD_IIb         = [];
                feature_FCD_IIb_cont    = [];
                
                for i = 1 : length(feature_FCD_IIb_idx)
                    
                    feature_FCD_IIb             = [ feature_FCD_IIb         mean(mean_z_lesion_fMRI_(feature_FCD_IIb_set{i}, FCD_IIb_fMRI), 1)' ];
                    if(length(feature_FCD_IIb_set{i}) > 1)
                        feature_FCD_IIb_cont    = [ feature_FCD_IIb_cont    mean(squeeze(mean(mean_z_lesion_fMRI_cont_(feature_FCD_IIb_set{i}, FCD_IIb_fMRI, :), 2)), 1)' ];
                    else
                        feature_FCD_IIb_cont    = [ feature_FCD_IIb_cont    mean(squeeze(mean_z_lesion_fMRI_cont_(feature_FCD_IIb_set{i}, FCD_IIb_fMRI, :)), 1)' ];
                    end
                    
                end
                
                abnorm = sum(sum(abs(feature_FCD_IIb)>IIb_threshold, 2)>=1);
                disp(['n=' num2str(abnorm), ' (' num2str(abnorm/sum(FCD_IIb_fMRI)) '%)']);
                
                %% multivariate features
                mahal_IIb      = mahal(feature_FCD_IIb, feature_FCD_IIb_cont);
                threshold_IIb  = mahal(IIb_threshold_set, feature_FCD_IIb_cont);
                disp(['n=' num2str(sum(mahal_IIb>=threshold_IIb)), ' (' num2str(sum(mahal_IIb>=threshold_IIb)/sum(FCD_IIb_fMRI)) '%)']);
                
                [a b] = max(mahal_IIb);
                case_num_pat_fMRI(FCDIIb_idx(b)) % highest md score: 065 071
                
                mahal_cont_IIb = [];
                for i = 1 : length(feature_FCD_IIb_cont)
                    feature_FCD_IIb_cont_temp = feature_FCD_IIb_cont;
                    feature_FCD_IIb_cont_temp(i, :) = [];
                    mahal_cont_IIb(i, :) = mahal(feature_FCD_IIb_cont(i, :), feature_FCD_IIb_cont_temp);
                end
                sum(mahal_cont_IIb>=threshold_IIb)/length(feature_FCD_IIb_cont)
                
            end
            
            %% All
            for All = 1
                threshold_sign = [ 1 -1 -1 -1 1 1 1 1 -1 -1 1 -1 1 -1 1 -1 1];
                IIb_threshold = 1.5;
                
                feature_FCD_IIb_idx = [2 3 4 5 7 8 9 11 15];
                
                feature_FCD_IIb_set = feature_set(feature_FCD_IIb_idx);
                IIb_threshold_sign = threshold_sign(feature_FCD_IIb_idx);
                IIb_threshold_set = IIb_threshold_sign*IIb_threshold;
                feature_FCD_IIb         = [];
                feature_FCD_IIb_cont    = [];
                
                for i = 1 : length(feature_FCD_IIb_idx)
                    
                    feature_FCD_IIb             = [ feature_FCD_IIb         mean(mean_z_lesion_fMRI_(feature_FCD_IIb_set{i}, FCD_IIb_fMRI), 1)' ];
                    if(length(feature_FCD_IIb_set{i}) > 1)
                        feature_FCD_IIb_cont    = [ feature_FCD_IIb_cont    mean(squeeze(mean(mean_z_lesion_fMRI_cont_(feature_FCD_IIb_set{i}, FCD_IIb_fMRI, :), 2)), 1)' ];
                    else
                        feature_FCD_IIb_cont    = [ feature_FCD_IIb_cont    mean(squeeze(mean_z_lesion_fMRI_cont_(feature_FCD_IIb_set{i}, FCD_IIb_fMRI, :)), 1)' ];
                    end
                    
                end
                
                mahal_IIb      = mahal(feature_FCD_IIb, feature_FCD_IIb_cont);
                threshold_IIb  = mahal(IIb_threshold_set, feature_FCD_IIb_cont);
                disp(['n=' num2str(sum(mahal_IIb>=threshold_IIb)), ' (' num2str(sum(mahal_IIb>=threshold_IIb)/sum(FCD_IIb_fMRI)) '%)']);
                
                [a b] = max(mahal_IIb);
                case_num_pat_fMRI(FCDIIb_fMRI_idx(b)) % highest md score: 065 071
                
                mahal_cont_IIb = [];
                for i = 1 : length(feature_FCD_IIb_cont)
                    feature_FCD_IIb_cont_temp = feature_FCD_IIb_cont;
                    feature_FCD_IIb_cont_temp(i, :) = [];
                    mahal_cont_IIb(i, :) = mahal(feature_FCD_IIb_cont(i, :), feature_FCD_IIb_cont_temp);
                end
                sum(mahal_cont_IIb>=threshold_IIb)/length(feature_FCD_IIb_cont)
                
            end
            
        end
    end
    
    for by_features = 1
                          %   1      2      3      4       5      6      7  8   9     10      11      12     13      14       15      16       17       18        19    20  21
        feature_set     = { [2:5], [7:9]  [16], [12:14], [15], [20:22], 17, 18, 19, [30:33], [35:37], [44], [40:42], [43], [45:47], [56 58], [60 62], [68 70], [72 74], 78, 79 };
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % IIa_individual_analysis
        for IIa_individual_analysis = 1
            
            threshold_sign = [ 1 -1 -1 -1 -1  -1 1 1 1 1 -1 -1 -1 -1 1 -1 -1 1 1 -1 1];
            
            %% Morphology
            for Morphology = 1
                
                IIa_threshold = 1.45;
                
                % each feature
                feature_FCD_IIa_idx = [ 7 8 9 ];
                
                feature_FCD_IIa_set = feature_set(feature_FCD_IIa_idx);
                IIa_threshold_sign = threshold_sign(feature_FCD_IIa_idx);
                IIa_threshold_set = IIa_threshold_sign*IIa_threshold;
                feature_FCD_IIa         = [];
                feature_FCD_IIa_cont    = [];
                
                for i = 1 : length(feature_FCD_IIa_idx)
                    
                    feature_FCD_IIa             = [ feature_FCD_IIa         mean(mean_z_lesion_(feature_FCD_IIa_set{i}, FCD_IIa_T1_FLAIR), 1)' ];
                    if(length(feature_FCD_IIa_set{i}) > 1)
                        feature_FCD_IIa_cont    = [ feature_FCD_IIa_cont    mean(squeeze(mean(mean_z_lesion_cont_(feature_FCD_IIa_set{i}, FCD_IIa_T1_FLAIR, :), 2)), 1)' ];
                    else
                        feature_FCD_IIa_cont    = [ feature_FCD_IIa_cont    mean(squeeze(mean_z_lesion_cont_(feature_FCD_IIa_set{i}, FCD_IIa_T1_FLAIR, :)), 1)' ];
                    end
                    
                end
                
                abnorm = sum(sum(abs(feature_FCD_IIa)>IIa_threshold, 2)>=1);
                disp(['n=' num2str(abnorm), ' (' num2str(abnorm/sum(FCD_IIa_T1_FLAIR)) '%)']);
                
                abnorm = sum(sum(abs(feature_FCD_IIa_cont)>IIa_threshold, 2)>=1);
                disp(['n=' num2str(abnorm), ' (' num2str(abnorm/length(feature_FCD_IIa_cont)) '%)']);

            end
            
            %% Intensity
            for Intensity = 1
                
                IIa_threshold = 1.45;
                
                % each feature
                feature_FCD_IIa_idx = [ 1 10 ]; % 1.5
                feature_FCD_IIa_idx = [ 6 15 ]; % 1.6
                feature_FCD_IIa_idx = [ 2 11 ]; % 1.5
                feature_FCD_IIa_idx = [ 3 12 ]; % 1.5
                feature_FCD_IIa_idx = [ 4 13 ]; % 1.5
                feature_FCD_IIa_idx = [ 5 14 ];
                
                feature_FCD_IIa_set = feature_set(feature_FCD_IIa_idx);
                IIa_threshold_sign = threshold_sign(feature_FCD_IIa_idx);
                IIa_threshold_set = IIa_threshold_sign*IIa_threshold;
                feature_FCD_IIa         = [];
                feature_FCD_IIa_cont    = [];
                
                for i = 1 : length(feature_FCD_IIa_idx)
                    
                    feature_FCD_IIa             = [ feature_FCD_IIa         mean(mean_z_lesion_(feature_FCD_IIa_set{i}, FCD_IIa_T1_FLAIR), 1)' ];
                    if(length(feature_FCD_IIa_set{i}) > 1)
                        feature_FCD_IIa_cont    = [ feature_FCD_IIa_cont    mean(squeeze(mean(mean_z_lesion_cont_(feature_FCD_IIa_set{i}, FCD_IIa_T1_FLAIR, :), 2)), 1)' ];
                    else
                        feature_FCD_IIa_cont    = [ feature_FCD_IIa_cont    mean(squeeze(mean_z_lesion_cont_(feature_FCD_IIa_set{i}, FCD_IIa_T1_FLAIR, :)), 1)' ];
                    end
                    
                end
                
                abnorm = sum(sum(abs(feature_FCD_IIa)>IIa_threshold, 2)>=1);
                disp(['n=' num2str(abnorm), ' (' num2str(abnorm/sum(FCD_IIa_T1_FLAIR)) '%)']);        
                
                abnorm = sum(sum(abs(feature_FCD_IIa_cont)>IIa_threshold, 2)>=1);
                disp(['n=' num2str(abnorm), ' (' num2str(abnorm/length(feature_FCD_IIa_cont)) '%)']);
               
            end
            
            %% DTI
            for DTI = 1
                
                IIa_threshold = 1.45;
                
                % each feature
                feature_FCD_IIa_idx = [ 16 18 ];
                feature_FCD_IIa_idx = [ 17 19 ];
                
                feature_FCD_IIa_set = feature_set(feature_FCD_IIa_idx);
                IIa_threshold_sign = threshold_sign(feature_FCD_IIa_idx);
                IIa_threshold_set = IIa_threshold_sign*IIa_threshold;
                feature_FCD_IIa         = [];
                feature_FCD_IIa_cont    = [];
                
                for i = 1 : length(feature_FCD_IIa_idx)
                    
                    feature_FCD_IIa             = [ feature_FCD_IIa         mean(mean_z_lesion_DTI_(feature_FCD_IIa_set{i}, FCD_IIa_DTI), 1)' ];
                    if(length(feature_FCD_IIa_set{i}) > 1)
                        feature_FCD_IIa_cont    = [ feature_FCD_IIa_cont    mean(squeeze(mean(mean_z_lesion_DTI_cont_(feature_FCD_IIa_set{i}, FCD_IIa_DTI, :), 2)), 1)' ];
                    else
                        feature_FCD_IIa_cont    = [ feature_FCD_IIa_cont    mean(squeeze(mean_z_lesion_DTI_cont_(feature_FCD_IIa_set{i}, FCD_IIa_DTI, :)), 1)' ];
                    end
                    
                end
                
                abnorm = sum(sum(abs(feature_FCD_IIa)>IIa_threshold, 2)>=1);
                disp(['n=' num2str(abnorm), ' (' num2str(abnorm/sum(FCD_IIa_DTI)) '%)']);        
                
                abnorm = sum(sum(abs(feature_FCD_IIa_cont)>IIa_threshold, 2)>=1);
                disp(['n=' num2str(abnorm), ' (' num2str(abnorm/length(feature_FCD_IIa_cont)) '%)']);
               
            end
            
            %% fMRI
            for fMRI = 1
               
                IIa_threshold = 1.45;
                
                % each feature
                feature_FCD_IIa_idx = [ 20 21 ];
                
                feature_FCD_IIa_set = feature_set(feature_FCD_IIa_idx);
                IIa_threshold_sign = threshold_sign(feature_FCD_IIa_idx);
                IIa_threshold_set = IIa_threshold_sign*IIa_threshold;
                feature_FCD_IIa         = [];
                feature_FCD_IIa_cont    = [];
                
                for i = 1 : length(feature_FCD_IIa_idx)
                    
                    feature_FCD_IIa             = [ feature_FCD_IIa         mean(mean_z_lesion_fMRI_(feature_FCD_IIa_set{i}, FCD_IIa_fMRI), 1)' ];
                    if(length(feature_FCD_IIa_set{i}) > 1)
                        feature_FCD_IIa_cont    = [ feature_FCD_IIa_cont    mean(squeeze(mean(mean_z_lesion_fMRI_cont_(feature_FCD_IIa_set{i}, FCD_IIa_fMRI, :), 2)), 1)' ];
                    else
                        feature_FCD_IIa_cont    = [ feature_FCD_IIa_cont    mean(squeeze(mean_z_lesion_fMRI_cont_(feature_FCD_IIa_set{i}, FCD_IIa_fMRI, :)), 1)' ];
                    end
                    
                end
                
                abnorm = sum(sum(abs(feature_FCD_IIa)>IIa_threshold, 2)>=1);
                disp(['n=' num2str(abnorm), ' (' num2str(abnorm/sum(FCD_IIa_fMRI)) '%)']);
                
                abnorm = sum(sum(abs(feature_FCD_IIa_cont)>IIa_threshold, 2)>=1);
                disp(['n=' num2str(abnorm), ' (' num2str(abnorm/length(feature_FCD_IIa_cont)) '%)']);
                
            
            end
            
            %% All
            for All = 1
                
                IIa_threshold = 1.45;
                
                feature_FCD_IIa_idx = [2 3 4 5 6 7 9 11 12 13 14 15 17 ];
                
                feature_FCD_IIa_set = feature_set(feature_FCD_IIa_idx);
                IIa_threshold_sign = threshold_sign(feature_FCD_IIa_idx);
                IIa_threshold_set = IIa_threshold_sign*IIa_threshold;
                feature_FCD_IIa         = [];
                feature_FCD_IIa_cont    = [];
                
                for i = 1 : length(feature_FCD_IIa_idx)
                    
                    feature_FCD_IIa             = [ feature_FCD_IIa         mean(mean_z_lesion_fMRI_(feature_FCD_IIa_set{i}, FCD_IIa_fMRI), 1)' ];
                    if(length(feature_FCD_IIa_set{i}) > 1)
                        feature_FCD_IIa_cont    = [ feature_FCD_IIa_cont    mean(squeeze(mean(mean_z_lesion_fMRI_cont_(feature_FCD_IIa_set{i}, FCD_IIa_fMRI, :), 2)), 1)' ];
                    else
                        feature_FCD_IIa_cont    = [ feature_FCD_IIa_cont    mean(squeeze(mean_z_lesion_fMRI_cont_(feature_FCD_IIa_set{i}, FCD_IIa_fMRI, :)), 1)' ];
                    end
                    
                end
                
                mahal_IIa      = mahal(feature_FCD_IIa, feature_FCD_IIa_cont);
                threshold_IIa  = mahal(IIa_threshold_set, feature_FCD_IIa_cont);
                disp(['n=' num2str(sum(mahal_IIa>=threshold_IIa)), ' (' num2str(sum(mahal_IIa>=threshold_IIa)/sum(FCD_IIa_fMRI)) '%)']);
                
                [a b] = max(mahal_IIa);
                case_num_pat_fMRI(FCDIIa_fMRI_idx(b)) % highest md score: 065 071
                
                mahal_cont_IIa = [];
                for i = 1 : length(feature_FCD_IIa_cont)
                    feature_FCD_IIa_cont_temp = feature_FCD_IIa_cont;
                    feature_FCD_IIa_cont_temp(i, :) = [];
                    mahal_cont_IIa(i, :) = mahal(feature_FCD_IIa_cont(i, :), feature_FCD_IIa_cont_temp);
                end
                sum(mahal_cont_IIa>=threshold_IIa)/length(feature_FCD_IIa_cont)
                
            end
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % IIb_individual_analysis
        for IIb_individual_analysis = 1
            
            threshold_sign = [ 1 -1 -1 -1 -1  -1 1 1 1 1 -1 -1 -1 -1 1 1 1 1 1 -1 -1];
            
            %% Morphology
            for Morphology = 1
                
                IIb_threshold = 1.45;
                
                % each feature
                feature_FCD_IIb_idx = [ 7 8 9 ];
                
                feature_FCD_IIb_set = feature_set(feature_FCD_IIb_idx);
                IIb_threshold_sign = threshold_sign(feature_FCD_IIb_idx);
                IIb_threshold_set = IIb_threshold_sign*IIb_threshold;
                feature_FCD_IIb         = [];
                feature_FCD_IIb_cont    = [];
                
                for i = 1 : length(feature_FCD_IIb_idx)
                    
                    feature_FCD_IIb             = [ feature_FCD_IIb         mean(mean_z_lesion_(feature_FCD_IIb_set{i}, FCD_IIb_T1_FLAIR), 1)' ];
                    if(length(feature_FCD_IIb_set{i}) > 1)
                        feature_FCD_IIb_cont    = [ feature_FCD_IIb_cont    mean(squeeze(mean(mean_z_lesion_cont_(feature_FCD_IIb_set{i}, FCD_IIb_T1_FLAIR, :), 2)), 1)' ];
                    else
                        feature_FCD_IIb_cont    = [ feature_FCD_IIb_cont    mean(squeeze(mean_z_lesion_cont_(feature_FCD_IIb_set{i}, FCD_IIb_T1_FLAIR, :)), 1)' ];
                    end
                    
                end
                
                abnorm = sum(sum(abs(feature_FCD_IIb)>IIb_threshold, 2)>=1);
                disp(['n=' num2str(abnorm), ' (' num2str(abnorm/sum(FCD_IIb_T1_FLAIR)) '%)']);
                
                abnorm = sum(sum(abs(feature_FCD_IIb_cont)>IIb_threshold, 2)>=1);
                disp(['n=' num2str(abnorm), ' (' num2str(abnorm/length(feature_FCD_IIb_cont)) '%)']);

            end
            
            %% Intensity
            for Intensity = 1
               
                IIb_threshold = 1.45;
                
                % each feature
                feature_FCD_IIb_idx = [ 1 10 ];
                feature_FCD_IIb_idx = [ 6 15 ];
                feature_FCD_IIb_idx = [ 2 11 ]; % 1.4
                feature_FCD_IIb_idx = [ 3 12 ];
                feature_FCD_IIb_idx = [ 4 13 ];
                feature_FCD_IIb_idx = [ 5 14 ];
                
                feature_FCD_IIb_set = feature_set(feature_FCD_IIb_idx);
                IIb_threshold_sign = threshold_sign(feature_FCD_IIb_idx);
                IIb_threshold_set = IIb_threshold_sign*IIb_threshold;
                feature_FCD_IIb         = [];
                feature_FCD_IIb_cont    = [];
                
                for i = 1 : length(feature_FCD_IIb_idx)
                    
                    feature_FCD_IIb             = [ feature_FCD_IIb         mean(mean_z_lesion_(feature_FCD_IIb_set{i}, FCD_IIb_T1_FLAIR), 1)' ];
                    if(length(feature_FCD_IIb_set{i}) > 1)
                        feature_FCD_IIb_cont    = [ feature_FCD_IIb_cont    mean(squeeze(mean(mean_z_lesion_cont_(feature_FCD_IIb_set{i}, FCD_IIb_T1_FLAIR, :), 2)), 1)' ];
                    else
                        feature_FCD_IIb_cont    = [ feature_FCD_IIb_cont    mean(squeeze(mean_z_lesion_cont_(feature_FCD_IIb_set{i}, FCD_IIb_T1_FLAIR, :)), 1)' ];
                    end
                    
                end
                
                abnorm = sum(sum(abs(feature_FCD_IIb)>IIb_threshold, 2)>=1);
                disp(['n=' num2str(abnorm), ' (' num2str(abnorm/sum(FCD_IIb_T1_FLAIR)) '%)']);        
                
                abnorm = sum(sum(abs(feature_FCD_IIb_cont)>IIb_threshold, 2)>=1);
                disp(['n=' num2str(abnorm), ' (' num2str(abnorm/length(feature_FCD_IIb_cont)) '%)']);
               
            end
            
            %% DTI
            for DTI = 1
                
                IIb_threshold = 1.45;
                
                % each feature
                feature_FCD_IIb_idx = [ 16 18 ];
                feature_FCD_IIb_idx = [ 17 19 ];
                
                feature_FCD_IIb_set = feature_set(feature_FCD_IIb_idx);
                IIb_threshold_sign = threshold_sign(feature_FCD_IIb_idx);
                IIb_threshold_set = IIb_threshold_sign*IIb_threshold;
                feature_FCD_IIb         = [];
                feature_FCD_IIb_cont    = [];
                
                for i = 1 : length(feature_FCD_IIb_idx)
                    
                    feature_FCD_IIb             = [ feature_FCD_IIb         mean(mean_z_lesion_DTI_(feature_FCD_IIb_set{i}, FCD_IIb_DTI), 1)' ];
                    if(length(feature_FCD_IIb_set{i}) > 1)
                        feature_FCD_IIb_cont    = [ feature_FCD_IIb_cont    mean(squeeze(mean(mean_z_lesion_DTI_cont_(feature_FCD_IIb_set{i}, FCD_IIb_DTI, :), 2)), 1)' ];
                    else
                        feature_FCD_IIb_cont    = [ feature_FCD_IIb_cont    mean(squeeze(mean_z_lesion_DTI_cont_(feature_FCD_IIb_set{i}, FCD_IIb_DTI, :)), 1)' ];
                    end
                    
                end
                
                abnorm = sum(sum(abs(feature_FCD_IIb)>IIb_threshold, 2)>=1);
                disp(['n=' num2str(abnorm), ' (' num2str(abnorm/sum(FCD_IIb_DTI)) '%)']);        
                
                abnorm = sum(sum(abs(feature_FCD_IIb_cont)>IIb_threshold, 2)>=1);
                disp(['n=' num2str(abnorm), ' (' num2str(abnorm/length(feature_FCD_IIb_cont)) '%)']);
               
            end
            
            %% fMRI
            for fMRI = 1
                
                IIb_threshold = 1.45;
                
                % each feature
                feature_FCD_IIb_idx = [ 20 21 ];
                
                feature_FCD_IIb_set = feature_set(feature_FCD_IIb_idx);
                IIb_threshold_sign = threshold_sign(feature_FCD_IIb_idx);
                IIb_threshold_set = IIb_threshold_sign*IIb_threshold;
                feature_FCD_IIb         = [];
                feature_FCD_IIb_cont    = [];
                
                for i = 1 : length(feature_FCD_IIb_idx)
                    
                    feature_FCD_IIb             = [ feature_FCD_IIb         mean(mean_z_lesion_fMRI_(feature_FCD_IIb_set{i}, FCD_IIb_fMRI), 1)' ];
                    if(length(feature_FCD_IIb_set{i}) > 1)
                        feature_FCD_IIb_cont    = [ feature_FCD_IIb_cont    mean(squeeze(mean(mean_z_lesion_fMRI_cont_(feature_FCD_IIb_set{i}, FCD_IIb_fMRI, :), 2)), 1)' ];
                    else
                        feature_FCD_IIb_cont    = [ feature_FCD_IIb_cont    mean(squeeze(mean_z_lesion_fMRI_cont_(feature_FCD_IIb_set{i}, FCD_IIb_fMRI, :)), 1)' ];
                    end
                    
                end
                
                abnorm = sum(sum(abs(feature_FCD_IIb)>IIb_threshold, 2)>=1);
                disp(['n=' num2str(abnorm), ' (' num2str(abnorm/sum(FCD_IIb_fMRI)) '%)']);
                
                abnorm = sum(sum(abs(feature_FCD_IIa_cont)>IIa_threshold, 2)>=1);
                disp(['n=' num2str(abnorm), ' (' num2str(abnorm/length(feature_FCD_IIa_cont)) '%)']);
                
            end
            
            %% All
            for All = 1
                
                IIb_threshold = 1.45;
                
                feature_FCD_IIb_idx = [2 3 5 6 7 9 11 12 15 19 20 21];
                
                feature_FCD_IIb_set = feature_set(feature_FCD_IIb_idx);
                IIb_threshold_sign = threshold_sign(feature_FCD_IIb_idx);
                IIb_threshold_set = IIb_threshold_sign*IIb_threshold;
                feature_FCD_IIb         = [];
                feature_FCD_IIb_cont    = [];
                
                for i = 1 : length(feature_FCD_IIb_idx)
                    
                    feature_FCD_IIb             = [ feature_FCD_IIb         mean(mean_z_lesion_fMRI_(feature_FCD_IIb_set{i}, FCD_IIb_fMRI), 1)' ];
                    if(length(feature_FCD_IIb_set{i}) > 1)
                        feature_FCD_IIb_cont    = [ feature_FCD_IIb_cont    mean(squeeze(mean(mean_z_lesion_fMRI_cont_(feature_FCD_IIb_set{i}, FCD_IIb_fMRI, :), 2)), 1)' ];
                    else
                        feature_FCD_IIb_cont    = [ feature_FCD_IIb_cont    mean(squeeze(mean_z_lesion_fMRI_cont_(feature_FCD_IIb_set{i}, FCD_IIb_fMRI, :)), 1)' ];
                    end
                    
                end
                
                mahal_IIb      = mahal(feature_FCD_IIb, feature_FCD_IIb_cont);
                threshold_IIb  = mahal(IIb_threshold_set, feature_FCD_IIb_cont);
                disp(['n=' num2str(sum(mahal_IIb>=threshold_IIb)), ' (' num2str(sum(mahal_IIb>=threshold_IIb)/sum(FCD_IIb_fMRI)) '%)']);
                
                [a b] = max(mahal_IIb);
                case_num_pat_fMRI(FCDIIb_fMRI_idx(b)) % highest md score: 065 071
                
                mahal_cont_IIb = [];
                for i = 1 : length(feature_FCD_IIb_cont)
                    feature_FCD_IIb_cont_temp = feature_FCD_IIb_cont;
                    feature_FCD_IIb_cont_temp(i, :) = [];
                    mahal_cont_IIb(i, :) = mahal(feature_FCD_IIb_cont(i, :), feature_FCD_IIb_cont_temp);
                end
                sum(mahal_cont_IIb>=threshold_IIb)/length(feature_FCD_IIb_cont)
                
            end
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Healthy controls_individual_analysis
        for Controls_individual_analysis = 1
            
            threshold_sign = [ 1 -1 -1 -1 -1  -1 1 1 1 1 -1 -1 -1 -1 1 1 1 1 1 -1 -1];
            
            %% Morphology
            for Morphology = 1
                
                IIb_threshold = 1.45;
                
                % each feature
                feature_cont_idx = [ 7 8 9 ];
                
                feature_cont_set = feature_set(feature_cont_idx);
                IIb_threshold_sign = threshold_sign(feature_cont_idx);
                IIb_threshold_set = IIb_threshold_sign*IIb_threshold;
                feature_FCD_IIb         = [];
                feature_FCD_IIb_cont    = [];
                
                for i = 1 : length(feature_cont_idx)
                                        
                    if(length(feature_cont_set{i}) > 1)
                        feature_FCD_IIb_cont    = [ feature_FCD_IIb_cont    mean(squeeze(mean(mean_z_lesion_cont_(feature_cont_set{i}, :, :), 2)), 1)' ];
                    else
                        feature_FCD_IIb_cont    = [ feature_FCD_IIb_cont    mean(squeeze(mean_z_lesion_cont_(feature_cont_set{i}, :, :)), 1)' ];
                    end
                    
                end
                
                abnorm = sum(sum(abs(feature_FCD_IIb_cont)>IIb_threshold, 2)>=1);
                disp(['n=' num2str(abnorm), ' (' num2str(abnorm/length(feature_FCD_IIb_cont)) '%)']);

            end
            
            %% Intensity
            for Intensity = 1
               
                IIb_threshold = 1.45;
                
                % each feature
                feature_cont_idx = [ 1 10 ];
                feature_cont_idx = [ 6 15 ];
                feature_cont_idx = [ 2 11 ]; % 1.4
                feature_cont_idx = [ 3 12 ];
                feature_cont_idx = [ 4 13 ];
                feature_cont_idx = [ 5 14 ];
                
                feature_cont_set = feature_set(feature_cont_idx);
                IIb_threshold_sign = threshold_sign(feature_cont_idx);
                IIb_threshold_set = IIb_threshold_sign*IIb_threshold;
                feature_FCD_IIb         = [];
                feature_FCD_IIb_cont    = [];
                
                for i = 1 : length(feature_cont_idx)
                                        
                    if(length(feature_cont_set{i}) > 1)
                        feature_FCD_IIb_cont    = [ feature_FCD_IIb_cont    mean(squeeze(mean(mean_z_lesion_cont_(feature_cont_set{i}, :, :), 2)), 1)' ];
                    else
                        feature_FCD_IIb_cont    = [ feature_FCD_IIb_cont    mean(squeeze(mean_z_lesion_cont_(feature_cont_set{i}, :, :)), 1)' ];
                    end
                    
                end
                
                abnorm = sum(sum(abs(feature_FCD_IIb_cont)>IIb_threshold, 2)>=1);
                disp(['n=' num2str(abnorm), ' (' num2str(abnorm/length(feature_FCD_IIb_cont)) '%)']);
               
            end
            
            %% DTI
            for DTI = 1
                
                IIb_threshold = 1.45;
                
                % each feature
                feature_cont_idx = [ 16 18 ];
                feature_cont_idx = [ 17 19 ];
                
                feature_cont_set = feature_set(feature_cont_idx);
                IIb_threshold_sign = threshold_sign(feature_cont_idx);
                IIb_threshold_set = IIb_threshold_sign*IIb_threshold;
                feature_FCD_IIb         = [];
                feature_FCD_IIb_cont    = [];
                
                for i = 1 : length(feature_cont_idx)
                                        
                    if(length(feature_cont_set{i}) > 1)
                        feature_FCD_IIb_cont    = [ feature_FCD_IIb_cont    mean(squeeze(mean(mean_z_lesion_DTI_cont_(feature_cont_set{i}, FCD_IIb_DTI, :), 2)), 1)' ];
                    else
                        feature_FCD_IIb_cont    = [ feature_FCD_IIb_cont    mean(squeeze(mean_z_lesion_DTI_cont_(feature_cont_set{i}, FCD_IIb_DTI, :)), 1)' ];
                    end
                    
                end
      
                abnorm = sum(sum(abs(feature_FCD_IIb_cont)>IIb_threshold, 2)>=1);
                disp(['n=' num2str(abnorm), ' (' num2str(abnorm/length(feature_FCD_IIb_cont)) '%)']);
               
            end
            
            %% fMRI
            for fMRI = 1
                
                IIb_threshold = 1.45;
                
                % each feature
                feature_cont_idx = [ 20 21 ];
                
                feature_cont_set = feature_set(feature_cont_idx);
                IIb_threshold_sign = threshold_sign(feature_cont_idx);
                IIb_threshold_set = IIb_threshold_sign*IIb_threshold;
                feature_FCD_IIb         = [];
                feature_FCD_IIb_cont    = [];
                
                for i = 1 : length(feature_cont_idx)
                
                    if(length(feature_cont_set{i}) > 1)
                        feature_FCD_IIb_cont    = [ feature_FCD_IIb_cont    mean(squeeze(mean(mean_z_lesion_fMRI_cont_(feature_cont_set{i}, FCD_IIb_fMRI, :), 2)), 1)' ];
                    else
                        feature_FCD_IIb_cont    = [ feature_FCD_IIb_cont    mean(squeeze(mean_z_lesion_fMRI_cont_(feature_cont_set{i}, FCD_IIb_fMRI, :)), 1)' ];
                    end
                    
                end
              
                abnorm = sum(sum(abs(feature_FCD_IIb_cont)>IIb_threshold, 2)>=1);
                disp(['n=' num2str(abnorm), ' (' num2str(abnorm/length(feature_FCD_IIb_cont)) '%)']);
                
            end
            
            %% All
            for All = 1
                
                IIb_threshold = 1.45;
                
                feature_cont_idx = [2 3 5 6 7 9 11 12 15 19 20 21];
                
                feature_cont_set = feature_set(feature_cont_idx);
                IIb_threshold_sign = threshold_sign(feature_cont_idx);
                IIb_threshold_set = IIb_threshold_sign*IIb_threshold;
                feature_FCD_IIb         = [];
                feature_FCD_IIb_cont    = [];
                
                for i = 1 : length(feature_cont_idx)
                    
                    feature_FCD_IIb             = [ feature_FCD_IIb         mean(mean_z_lesion_fMRI_(feature_cont_set{i}, FCD_IIb_fMRI), 1)' ];
                    if(length(feature_cont_set{i}) > 1)
                        feature_FCD_IIb_cont    = [ feature_FCD_IIb_cont    mean(squeeze(mean(mean_z_lesion_fMRI_cont_(feature_cont_set{i}, FCD_IIb_fMRI, :), 2)), 1)' ];
                    else
                        feature_FCD_IIb_cont    = [ feature_FCD_IIb_cont    mean(squeeze(mean_z_lesion_fMRI_cont_(feature_cont_set{i}, FCD_IIb_fMRI, :)), 1)' ];
                    end
                    
                end
                
                mahal_IIb      = mahal(feature_FCD_IIb, feature_FCD_IIb_cont);
                threshold_IIb  = mahal(IIb_threshold_set, feature_FCD_IIb_cont);
              
                mahal_cont_IIb = [];
                for i = 1 : length(feature_FCD_IIb_cont)
                    feature_FCD_IIb_cont_temp = feature_FCD_IIb_cont;
                    feature_FCD_IIb_cont_temp(i, :) = [];
                    mahal_cont_IIb(i, :) = mahal(feature_FCD_IIb_cont(i, :), feature_FCD_IIb_cont_temp);
                end
                sum(mahal_cont_IIb>=threshold_IIb)/length(feature_FCD_IIb_cont)
                
            end
            
        end        
        
    end
    
end