clear all
close all
addpath('/data/noel/noel2/local/matlab/surfstat_chicago');
addpath('/local_raid/seokjun/03_downloads/LSSVMlabv1_8_R2009b_R2011a');

%% Read pre-calculated zscore of T1 and FLAIR in the lesion
for read_matfiles = 1
    
    %% global-depth lesion profiling
    load('zscore_database_T1_FLAIR.mat');
    mean_z_lesion_T1_FLAIR = mean_z_lesion;
    load('zscore_database_DTI.mat');
    mean_z_lesion_DTI = mean_z_lesion;
    load('zscore_database_fMRI.mat');
    mean_z_lesion_fMRI = mean_z_lesion;
    
    %% perilesional profiling
    Kernel = 2;    
    load(['pl_feature_set_sm_' num2str(Kernel) '_T1_FLAIR.mat']);
    pl_feature_set_T1_FLAIR = pl_feature_set;
    load(['pl_feature_set_sm_' num2str(Kernel) '_DTI.mat']);
    pl_feature_set_DTI = pl_feature_set;
    load(['pl_feature_set_sm_' num2str(Kernel) '_fMRI.mat']);
    pl_feature_set_fMRI = pl_feature_set;
    load('zscore_database_perilesional_T1_FLAIR_DTI_fMRI_prediction.mat');
    
    %% statistical test permutation result
    load('permutation_test_classifier.mat');
    
end

%% Read demograpic data and set up some parameters and variables
for read_demodata_setup_params = 1
    
    OUTPATH = '/local_raid/seokjun/01_project/04_IntracorticalAnalysis/03_result/prediction_histology_by_MRI/';
    Cases_cont                = '/local_raid/seokjun/01_project/04_IntracorticalAnalysis/01_analysis/Demographic_data_control.txt';
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
    
    fid = fopen(Cases_cont);
    demo = textscan(fid, '%s%f%s', 'Delimiter', ',', 'CollectOutput', 1);
    case_num_cont = demo{1};
    age_cont = demo{2};
    gender_cont = demo{3};
    fclose(fid);
    
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

%% Predict histological subtype based different combination of features
for t1_flair_dti_fMRI = 1
    
    %% Make the feature vectors, matching the T1, FLAIR and DTI zscore
    for featureset_build = 1
        
        [C, ia1, ib1] = intersect(case_num_pat, case_num_pat_fMRI);
        [C, ia2, ib2] = intersect(case_num_pat_dti, case_num_pat_fMRI);
        FCD_IIa_org = strncmp(histo_type, 'FCDIIa', 6)';
        FCD_IIb_org = strncmp(histo_type, 'FCDIIb', 6)';
        FCD_IIa_dti = strncmp(histo_type_dti, 'FCDIIa', 6)';
        FCD_IIb_dti = strncmp(histo_type_dti, 'FCDIIb', 6)';
        FCD_IIa_fMRI = strncmp(histo_type_fMRI, 'FCDIIa', 6)';
        FCD_IIb_fMRI = strncmp(histo_type_fMRI, 'FCDIIb', 6)';
        
        mean_z_perilesion_T1_final    = zeros(1, length(case_num_pat));
        mean_z_perilesion_FLAIR_final = zeros(6, length(case_num_pat));
        mean_z_perilesion_DTI_final   = zeros(5, length(case_num_pat_dti));
        mean_z_perilesion_fMRI_final  = zeros(14, length(case_num_pat_fMRI));
        
        mean_z_perilesion_T1_final(FCD_IIa_org) = mean_z_perilesion_T1.mean_z_score_typeIIa';
        mean_z_perilesion_T1_final(FCD_IIb_org) = mean_z_perilesion_T1.mean_z_score_typeIIb';
        mean_z_perilesion_FLAIR_final(:, FCD_IIa_org) = mean_z_perilesion_FLAIR.mean_z_score_typeIIa';
        mean_z_perilesion_FLAIR_final(:, FCD_IIb_org) = mean_z_perilesion_FLAIR.mean_z_score_typeIIb';
        mean_z_perilesion_DTI_final(:, FCD_IIa_dti) = mean_z_perilesion_DTI.mean_z_score_typeIIa';
        mean_z_perilesion_DTI_final(:, FCD_IIb_dti) = mean_z_perilesion_DTI.mean_z_score_typeIIb';
        mean_z_perilesion_fMRI_final(:, FCD_IIa_fMRI) = mean_z_perilesion_fMRI.mean_z_score_typeIIa';
        mean_z_perilesion_fMRI_final(:, FCD_IIb_fMRI) = mean_z_perilesion_fMRI.mean_z_score_typeIIb';
        
        mean_z_perilesion_T1_final = mean_z_perilesion_T1_final(1, ia1);
        mean_z_perilesion_FLAIR_final = mean_z_perilesion_FLAIR_final(:, ia1);
        mean_z_perilesion_DTI_final = mean_z_perilesion_DTI_final(:, ia2); 
        
        mean_z_lesion      = [ mean_z_lesion_T1_FLAIR(:, ia1); mean_z_lesion_DTI(:, ia2); mean_z_lesion_fMRI; mean(mean_z_perilesion_T1_final, 1); mean(mean_z_perilesion_FLAIR_final, 1); mean(mean_z_perilesion_DTI_final, 1); mean(mean_z_perilesion_fMRI_final, 1) ];
        mean_z_lesion_temp = [ mean_z_lesion_T1_FLAIR(:, ia1); mean_z_lesion_DTI(:, ia2); mean_z_lesion_fMRI; ];
        
        peri_lesional_dist      =  [ 2:2:4 ];
        peri_lesional_dist_fMRI =  [ 2:4:8 ];
        mean_z_lesion_temp2 = [ mean_z_lesion_T1_FLAIR(:, ia1); mean_z_lesion_DTI(:, ia2); mean_z_lesion_fMRI;
            pl_feature_set_T1_FLAIR(6).data(peri_lesional_dist, ia1);                                                                                          %% T1 cortical thickness
            mean(pl_feature_set_T1_FLAIR(9).data(peri_lesional_dist, ia1, 2:5), 3);                                                                            %% FLAIR cortical intentisy
            mean(pl_feature_set_DTI(1).data(peri_lesional_dist, ia2, [3 5]), 3);                                                                               %% DTI cortical FA
            pl_feature_set_fMRI(1).data(peri_lesional_dist_fMRI, :);                                                                                           %% fMRI ALFF
            pl_feature_set_fMRI(2).data(peri_lesional_dist_fMRI, :);                                                                                           %% fMRI ReHo
            ];
        
    end
    
    %% 1. Use global features based on pre-selection
    for prediction_using_unimodal_global_features = 1
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Feature set index
        %% mean_z_lesion
        %% 1:5   = RI of T1                      -> 2:5      p1
        %% 6:10  = PG of T1                      -> 7:9 16   p2
        %% 11:15 = TG of T1                      -> 12:15    p3
        %% 16    = PG_gw of T1                   -> 7:9 16   p2
        %% 17    = CT of T1                      -> 17       p7
        %% 18    = MC of T1                      -> 18       p8
        %% 19    = SD of T1                      -> 19       p9
        %% 20:22 = RI of T1 subcortical          -> 20:22    p4
        %% 23:25 = PG of T1 subcortical          -> 23:25    p5
        %% 26:28 = TG of T1 subcortical          -> 26:28    p6
        %% 29:33 = RI of FLAIR                   -> 30:33    p10
        %% 34:38 = PG of FLAIR                   -> 35:37 44 p11
        %% 39:43 = TG of FLAIR                   -> 40:43    p12
        %% 44    = PG_gw of FLAIR                -> 35:37 44 p11
        %% 45:47 = RI of FLAIR subcortical       -> 45:47    p13
        %% 48:50 = PG of FLAIR subcortical       -> 48:50    p14
        %% 51:53 = TG of FLAIR subcortical       -> 51:53    p15
        %% 54:58 = FA of DTI Intra               -> 56 58    p16
        %% 59:65 = FA of DT1 Subcortical         -> 60 62    p17
        %% 66:70 = MD of DTI Intra               -> 68 70    p18
        %% 71:77 = MD of DTI Subcortical         -> 72 74    p19
        %% 78    = ALFF of functional            -> 78       p20
        %% 79    = ReHo of functional            -> 79       p21
        %% 80    = DC of functional              -> 80       p22
        %% 81    = BC of functional              -> 81       p23
        %% 82    = Thickness of T1 perilesion    -> 82       p24
        %% 83    = Intensity of FLAIR perilesion -> 83       p25
        %% 84    = FA of DTI perilesion          -> 84       p26
        %% 85    = FA of fMRI perilesion         -> 85       p27
        
        % Using t1 global features
        for T1_feature_gb_only = 1
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Feature set index
            %% mean_z_lesion
            %% 1:5   = RI of T1                      -> 2:5      p1
            %% 6:10  = PG of T1                      -> 7:9 16   p2
            %% 11:15 = TG of T1                      -> 12:15    p3
            %% 16    = PG_gw of T1                   -> 7:9 16   p2
            %% 17    = CT of T1                      -> 17       p7
            %% 18    = MC of T1                      -> 18       p8
            %% 19    = SD of T1                      -> 19       p9
            %% 20:22 = RI of T1 subcortical          -> 20:22    p4
            %% 23:25 = PG of T1 subcortical          -> 23:25    p5
            %% 26:28 = TG of T1 subcortical          -> 26:28    p6
            
            discriminative_features = { [2:5], [7:9 16], [12:15], 17, 18, 19, [20:22] };
            mean_z_lesion_final = mean_z_lesion;
            classified_histo = cell(length(histo_type_fMRI), 1);
            SVMstruct_set = cell(length(histo_type_fMRI), 1);
            discriminative_features_set = [];
            
            opts = statset('display', 'iter');            
           
            parfor i = 1 : length(histo_type_fMRI)
                
                LOOCV_sub = 1 : length(histo_type_fMRI);
                LOOCV_sub(i) = [];
                LOOCV_group = histo_type_fMRI;
                LOOCV_group(i) = [];
                                
                mean_z_lesion_temp_lvo = mean_z_lesion_final;                
                mean_z_lesion_temp_lvo2 = [];
                mean_z_lesion_temp_lvo3 = [];
                
                for f = 1 : length(discriminative_features)
                    mean_z_lesion_temp_lvo2 = [ mean_z_lesion_temp_lvo2; mean(mean_z_lesion_temp_lvo(discriminative_features{f}, LOOCV_sub), 1) ];
                    mean_z_lesion_temp_lvo3 = [ mean_z_lesion_temp_lvo3; mean(mean_z_lesion_temp_lvo(discriminative_features{f}, i), 1) ];
                end
                
                c = cvpartition(LOOCV_group,'resubstitution');
                [fs,history] = sequentialfs(@fs_func_svm,mean_z_lesion_temp_lvo2',LOOCV_group,'cv',c, 'direction', 'forward', 'options',opts);
                
                discriminative_features_fs = fs;
                discriminative_features_set = [ discriminative_features_set; discriminative_features_fs ];
                training_data = mean_z_lesion_temp_lvo2(discriminative_features_fs, :);
                new_data = mean_z_lesion_temp_lvo3(discriminative_features_fs); 
                
                SVMstruct = svmtrain(training_data', LOOCV_group, 'Kernel_Function', 'polynomial');
                newClass = svmclassify(SVMstruct, new_data');
                classified_histo(i) = newClass;
                SVMstruct_set(i) = { SVMstruct };
                
            end
            predict_T1w_global = (sum(strcmp(classified_histo, histo_type_fMRI)))/length(histo_type_fMRI); 
            
            %% A diagnostic test
            diagnostic_test_T1w_global = [
                sum(strcmp(histo_type_fMRI(strcmp(histo_type_fMRI, 'FCDIIa'), 1), classified_histo(strcmp(histo_type_fMRI, 'FCDIIa'), 1))) sum(~(strcmp(histo_type_fMRI(strcmp(histo_type_fMRI, 'FCDIIb'), 1), classified_histo(strcmp(histo_type_fMRI, 'FCDIIb'), 1))));
                sum(~(strcmp(histo_type_fMRI(strcmp(histo_type_fMRI, 'FCDIIa'), 1), classified_histo(strcmp(histo_type_fMRI, 'FCDIIa'), 1)))) sum(strcmp(histo_type_fMRI(strcmp(histo_type_fMRI, 'FCDIIb'), 1), classified_histo(strcmp(histo_type_fMRI, 'FCDIIb'), 1)))
                ];        
            
            discriminative_features_set_T1_SingleGlobal = discriminative_features_set;
            sum(discriminative_features_set, 1) % Vertical gradient and sulcal depth
            predict_T1w_global_case = strcmp(classified_histo, histo_type_fMRI);
            
            %% ------------------------------------------------
            %%                               Histology
            %%                       FCD IIa         FCD IIb
            %% ------------------------------------------------
            %%            FCD IIa      5        |       8
            %% Prediction -------------------------------------
            %%            FCD IIb      3        |       7
            %% ------------------------------------------------
            %% sensitivity: 12/23 -> 52%
            %%        IIa : 5/8   -> 63%
            %%        IIb : 7/15  -> 47%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        end
        
        % Using FLAIR global features
        for FLAIR_feature_gb_only = 1
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Feature set index
            %% mean_z_lesion            
            %% 29:33 = RI of FLAIR                   -> 30:33    p10
            %% 34:38 = PG of FLAIR                   -> 35:37 44 p11
            %% 39:43 = TG of FLAIR                   -> 40:43    p12
            %% 44    = PG_gw of FLAIR                -> 35:37 44 p11
            %% 45:47 = RI of FLAIR subcortical       -> 45:47    p13
            %% 48:50 = PG of FLAIR subcortical       -> 48:50    p14
            %% 51:53 = TG of FLAIR subcortical       -> 51:53    p15
            
            discriminative_features = { [30:33], [35:37 44], [40:43], [45:47] };
            mean_z_lesion_final = mean_z_lesion;
            classified_histo = cell(length(histo_type_fMRI), 1);
            SVMstruct_set = cell(length(histo_type_fMRI), 1);
            discriminative_features_set = [];
            
            opts = statset('display', 'iter');            
           
            parfor i = 1 : length(histo_type_fMRI)
                
                LOOCV_sub = 1 : length(histo_type_fMRI);
                LOOCV_sub(i) = [];
                LOOCV_group = histo_type_fMRI;
                LOOCV_group(i) = [];
                                
                mean_z_lesion_temp_lvo = mean_z_lesion_final;                
                mean_z_lesion_temp_lvo2 = [];
                mean_z_lesion_temp_lvo3 = [];
                
                for f = 1 : length(discriminative_features)
                    mean_z_lesion_temp_lvo2 = [ mean_z_lesion_temp_lvo2; mean(mean_z_lesion_temp_lvo(discriminative_features{f}, LOOCV_sub), 1) ];
                    mean_z_lesion_temp_lvo3 = [ mean_z_lesion_temp_lvo3; mean(mean_z_lesion_temp_lvo(discriminative_features{f}, i), 1) ];
                end
                
                c = cvpartition(LOOCV_group,'resubstitution');
                [fs,history] = sequentialfs(@fs_func_svm,mean_z_lesion_temp_lvo2',LOOCV_group,'cv',c, 'direction', 'backward', 'options',opts);
                
                discriminative_features_fs = fs;
                discriminative_features_set = [ discriminative_features_set; discriminative_features_fs ];
                training_data = mean_z_lesion_temp_lvo2(discriminative_features_fs, :);
                new_data = mean_z_lesion_temp_lvo3(discriminative_features_fs); 
                
                SVMstruct = svmtrain(training_data', LOOCV_group, 'Kernel_Function', 'polynomial');
                newClass = svmclassify(SVMstruct, new_data');
                classified_histo(i) = newClass;
                SVMstruct_set(i) = { SVMstruct };
                
            end
            predict_FLAIR_global = (sum(strcmp(classified_histo, histo_type_fMRI)))/length(histo_type_fMRI);
            
            %% A diagnostic test
            diagnostic_test_FLAIR_global = [
                sum(strcmp(histo_type_fMRI(strcmp(histo_type_fMRI, 'FCDIIa'), 1), classified_histo(strcmp(histo_type_fMRI, 'FCDIIa'), 1))) sum(~(strcmp(histo_type_fMRI(strcmp(histo_type_fMRI, 'FCDIIb'), 1), classified_histo(strcmp(histo_type_fMRI, 'FCDIIb'), 1))));
                sum(~(strcmp(histo_type_fMRI(strcmp(histo_type_fMRI, 'FCDIIa'), 1), classified_histo(strcmp(histo_type_fMRI, 'FCDIIa'), 1)))) sum(strcmp(histo_type_fMRI(strcmp(histo_type_fMRI, 'FCDIIb'), 1), classified_histo(strcmp(histo_type_fMRI, 'FCDIIb'), 1)))
                ];
            
            discriminative_features_set_FLAIR_SingleGlobal = discriminative_features_set;
            sum(discriminative_features_set, 1) % Cortical intensity, vertical gradient and subcortical intensity
            predict_FLAIR_global_case = strcmp(classified_histo, histo_type_fMRI);
                        
            %% ------------------------------------------------
            %%                               Histology
            %%                       FCD IIa         FCD IIb
            %% ------------------------------------------------
            %%            FCD IIa      5        |       5
            %% Prediction -------------------------------------
            %%            FCD IIb      3        |       10
            %% ------------------------------------------------
            %% sensitivity: 15/23 -> 65%
            %%        IIa : 5/8   -> 63%
            %%        IIb : 10/15 -> 67%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        end
        
        % Using DTI global features
        for DTI_feature_gb_only = 1
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Feature set index
            %% mean_z_lesion
            %% 54:58 = FA of DTI Intra               -> 56 58    p16
            %% 59:65 = FA of DT1 Subcortical         -> 60 62    p17
            %% 66:70 = MD of DTI Intra               -> 68 70    p18
            %% 71:77 = MD of DTI Subcortical         -> 72 74    p19
            
            discriminative_features = { [56 58], [60 62], [68 70], [72 74] };
            mean_z_lesion_final = mean_z_lesion;
            classified_histo = cell(length(histo_type_fMRI), 1);
            SVMstruct_set = cell(length(histo_type_fMRI), 1);
            discriminative_features_set = [];
            
            opts = statset('display', 'iter');            
           
            parfor i = 1 : length(histo_type_fMRI)
                
                LOOCV_sub = 1 : length(histo_type_fMRI);
                LOOCV_sub(i) = [];
                LOOCV_group = histo_type_fMRI;
                LOOCV_group(i) = [];
                                
                mean_z_lesion_temp_lvo = mean_z_lesion_final;                
                mean_z_lesion_temp_lvo2 = [];
                mean_z_lesion_temp_lvo3 = [];
                
                for f = 1 : length(discriminative_features)
                    mean_z_lesion_temp_lvo2 = [ mean_z_lesion_temp_lvo2; mean(mean_z_lesion_temp_lvo(discriminative_features{f}, LOOCV_sub), 1) ];
                    mean_z_lesion_temp_lvo3 = [ mean_z_lesion_temp_lvo3; mean(mean_z_lesion_temp_lvo(discriminative_features{f}, i), 1) ];
                end
                
                c = cvpartition(LOOCV_group,'resubstitution');
                [fs,history] = sequentialfs(@fs_func_svm,mean_z_lesion_temp_lvo2',LOOCV_group,'cv',c, 'direction', 'forward', 'options',opts);
                
                discriminative_features_fs = fs;
                discriminative_features_set = [ discriminative_features_set; discriminative_features_fs ];
                training_data = mean_z_lesion_temp_lvo2(discriminative_features_fs, :);
                new_data = mean_z_lesion_temp_lvo3(discriminative_features_fs); 
                
                SVMstruct = svmtrain(training_data', LOOCV_group, 'Kernel_Function', 'polynomial');
                newClass = svmclassify(SVMstruct, new_data');
                classified_histo(i) = newClass;
                SVMstruct_set(i) = { SVMstruct };
                
            end
            predict_DTI_global = (sum(strcmp(classified_histo, histo_type_fMRI)))/length(histo_type_fMRI);
            
            %% A diagnostic test
            diagnostic_test_DTI_global = [
                sum(strcmp(histo_type_fMRI(strcmp(histo_type_fMRI, 'FCDIIa'), 1), classified_histo(strcmp(histo_type_fMRI, 'FCDIIa'), 1))) sum(~(strcmp(histo_type_fMRI(strcmp(histo_type_fMRI, 'FCDIIb'), 1), classified_histo(strcmp(histo_type_fMRI, 'FCDIIb'), 1))));
                sum(~(strcmp(histo_type_fMRI(strcmp(histo_type_fMRI, 'FCDIIa'), 1), classified_histo(strcmp(histo_type_fMRI, 'FCDIIa'), 1)))) sum(strcmp(histo_type_fMRI(strcmp(histo_type_fMRI, 'FCDIIb'), 1), classified_histo(strcmp(histo_type_fMRI, 'FCDIIb'), 1)))
                ];
            
            discriminative_features_set_DTI_SingleGlobal = discriminative_features_set;
            sum(discriminative_features_set, 1) % Cortical and subcortical FA
            predict_DTI_global_case = strcmp(classified_histo, histo_type_fMRI);
                        
            %% ------------------------------------------------
            %%                               Histology
            %%                       FCD IIa         FCD IIb
            %% ------------------------------------------------
            %%            FCD IIa      3        |       4
            %% Prediction -------------------------------------
            %%            FCD IIb      5        |       11
            %% ------------------------------------------------
            %% sensitivity: 14/23 -> 61%
            %%        IIa : 3/8   -> 38%
            %%        IIb : 11/15 -> 73%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        end
        
        % Using fMRI global features
        for fMRI_feature_gb_only = 1
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Feature set index
            %% mean_z_lesion
            %% 78    = ALFF of functional            -> 78       p20
            %% 79    = ReHo of functional            -> 79       p21
            %% 80    = DC of functional              -> 80       p22
            %% 81    = BC of functional              -> 81       p23
            
            discriminative_features = { 78, 79 };
            mean_z_lesion_final = mean_z_lesion;
            classified_histo = cell(length(histo_type_fMRI), 1);
            SVMstruct_set = cell(length(histo_type_fMRI), 1);
            discriminative_features_set = [];
            
            opts = statset('display', 'iter');            
           
            parfor i = 1 : length(histo_type_fMRI)
                
                LOOCV_sub = 1 : length(histo_type_fMRI);
                LOOCV_sub(i) = [];
                LOOCV_group = histo_type_fMRI;
                LOOCV_group(i) = [];
                                
                mean_z_lesion_temp_lvo = mean_z_lesion_final;                
                mean_z_lesion_temp_lvo2 = [];
                mean_z_lesion_temp_lvo3 = [];
                
                for f = 1 : length(discriminative_features)
                    mean_z_lesion_temp_lvo2 = [ mean_z_lesion_temp_lvo2; mean(mean_z_lesion_temp_lvo(discriminative_features{f}, LOOCV_sub), 1) ];
                    mean_z_lesion_temp_lvo3 = [ mean_z_lesion_temp_lvo3; mean(mean_z_lesion_temp_lvo(discriminative_features{f}, i), 1) ];
                end
                
                c = cvpartition(LOOCV_group,'resubstitution');
                [fs,history] = sequentialfs(@fs_func_svm,mean_z_lesion_temp_lvo2',LOOCV_group,'cv',c, 'direction', 'forward', 'options',opts);
                
                discriminative_features_fs = fs;
                discriminative_features_set = [ discriminative_features_set; discriminative_features_fs ];
                training_data = mean_z_lesion_temp_lvo2(discriminative_features_fs, :);
                new_data = mean_z_lesion_temp_lvo3(discriminative_features_fs); 
                
                SVMstruct = svmtrain(training_data', LOOCV_group, 'Kernel_Function', 'polynomial');
                newClass = svmclassify(SVMstruct, new_data');
                classified_histo(i) = newClass;
                SVMstruct_set(i) = { SVMstruct };
                
            end
            predict_fMRI_global=(sum(strcmp(classified_histo, histo_type_fMRI)))/length(histo_type_fMRI);
            
            %% A diagnostic test
            diagnostic_test_fMRI_global = [
                sum(strcmp(histo_type_fMRI(strcmp(histo_type_fMRI, 'FCDIIa'), 1), classified_histo(strcmp(histo_type_fMRI, 'FCDIIa'), 1))) sum(~(strcmp(histo_type_fMRI(strcmp(histo_type_fMRI, 'FCDIIb'), 1), classified_histo(strcmp(histo_type_fMRI, 'FCDIIb'), 1))));
                sum(~(strcmp(histo_type_fMRI(strcmp(histo_type_fMRI, 'FCDIIa'), 1), classified_histo(strcmp(histo_type_fMRI, 'FCDIIa'), 1)))) sum(strcmp(histo_type_fMRI(strcmp(histo_type_fMRI, 'FCDIIb'), 1), classified_histo(strcmp(histo_type_fMRI, 'FCDIIb'), 1)))
                ];
            
            discriminative_features_set_fMRI_SingleGlobal = discriminative_features_set;
            sum(discriminative_features_set, 1) % ALFF and ReHo
            predict_fMRI_global_case = strcmp(classified_histo, histo_type_fMRI);
                        
            %% ------------------------------------------------
            %%                               Histology
            %%                       FCD IIa         FCD IIb
            %% ------------------------------------------------
            %%            FCD IIa      8        |       5
            %% Prediction -------------------------------------
            %%            FCD IIb      1        |       9
            %% ------------------------------------------------
            %% sensitivity: 17/23 -> 74%
            %%        IIa : 8/9   -> 89%
            %%        IIb : 9/14  -> 64%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        end      
        
        % Using all global features
        for all_modalities_gb_only = 1
            
            discriminative_features = { [2:5], [7:9 16], [12:15], 17, 18, 19, [20:22], [30:33], [35:37 44], [40:43], [45:47], [56 58], [60 62], [68 70], [72 74], [78, 79] };
            mean_z_lesion_final = mean_z_lesion;
            classified_histo = cell(length(histo_type_fMRI), 1);
            SVMstruct_set = cell(length(histo_type_fMRI), 1);
            discriminative_features_set = [];
            
            opts = statset('display', 'iter');            
           
            parfor i = 1 : length(histo_type_fMRI)
                
                LOOCV_sub = 1 : length(histo_type_fMRI);
                LOOCV_sub(i) = [];
                LOOCV_group = histo_type_fMRI;
                LOOCV_group(i) = [];
                                
                mean_z_lesion_temp_lvo = mean_z_lesion_final;                
                mean_z_lesion_temp_lvo2 = [];
                mean_z_lesion_temp_lvo3 = [];
                
                for f = 1 : length(discriminative_features)
                    mean_z_lesion_temp_lvo2 = [ mean_z_lesion_temp_lvo2; mean(mean_z_lesion_temp_lvo(discriminative_features{f}, LOOCV_sub), 1) ];
                    mean_z_lesion_temp_lvo3 = [ mean_z_lesion_temp_lvo3; mean(mean_z_lesion_temp_lvo(discriminative_features{f}, i), 1) ];
                end
                
                c = cvpartition(LOOCV_group,'resubstitution');
                [fs,history] = sequentialfs(@fs_func_svm,mean_z_lesion_temp_lvo2',LOOCV_group,'cv',c, 'direction', 'forward', 'options',opts);
                
                discriminative_features_fs = fs;
                discriminative_features_set = [ discriminative_features_set; discriminative_features_fs ];
                training_data = mean_z_lesion_temp_lvo2(discriminative_features_fs, :);
                new_data = mean_z_lesion_temp_lvo3(discriminative_features_fs); 
                
                SVMstruct = svmtrain(training_data', LOOCV_group, 'Kernel_Function', 'polynomial');
                newClass = svmclassify(SVMstruct, new_data');
                classified_histo(i) = newClass;
                SVMstruct_set(i) = { SVMstruct };
                
            end
            predict_AllMod_global = (sum(strcmp(classified_histo, histo_type_fMRI)))/length(histo_type_fMRI);
            
            %% A diagnostic test
            diagnostic_test_AllMod_global = [
                sum(strcmp(histo_type_fMRI(strcmp(histo_type_fMRI, 'FCDIIa'), 1), classified_histo(strcmp(histo_type_fMRI, 'FCDIIa'), 1))) sum(~(strcmp(histo_type_fMRI(strcmp(histo_type_fMRI, 'FCDIIb'), 1), classified_histo(strcmp(histo_type_fMRI, 'FCDIIb'), 1))));
                sum(~(strcmp(histo_type_fMRI(strcmp(histo_type_fMRI, 'FCDIIa'), 1), classified_histo(strcmp(histo_type_fMRI, 'FCDIIa'), 1)))) sum(strcmp(histo_type_fMRI(strcmp(histo_type_fMRI, 'FCDIIb'), 1), classified_histo(strcmp(histo_type_fMRI, 'FCDIIb'), 1)))
                ];           
            
            discriminative_features_set_MultimodalGlobal = discriminative_features_set;
            sum(discriminative_features_set, 1) % Vertical gradient and ALFF + ReHo
            predict_AllMod_global_case = strcmp(classified_histo, histo_type_fMRI);
            
            % only structure global
            discriminative_features = { [2:5], [7:9 16], [12:15], 17, 18, 19, [20:22], [30:33], [35:37 44], [40:43], [45:47], [56 58], [60 62], [68 70], [72 74] };
            mean_z_lesion_final = mean_z_lesion;
            classified_histo = cell(length(histo_type_fMRI), 1);
            SVMstruct_set = cell(length(histo_type_fMRI), 1);
            discriminative_features_set = [];
            
            opts = statset('display', 'iter');            
           
            parfor i = 1 : length(histo_type_fMRI)
                
                LOOCV_sub = 1 : length(histo_type_fMRI);
                LOOCV_sub(i) = [];
                LOOCV_group = histo_type_fMRI;
                LOOCV_group(i) = [];
                                
                mean_z_lesion_temp_lvo = mean_z_lesion_final;                
                mean_z_lesion_temp_lvo2 = [];
                mean_z_lesion_temp_lvo3 = [];
                
                for f = 1 : length(discriminative_features)
                    mean_z_lesion_temp_lvo2 = [ mean_z_lesion_temp_lvo2; mean(mean_z_lesion_temp_lvo(discriminative_features{f}, LOOCV_sub), 1) ];
                    mean_z_lesion_temp_lvo3 = [ mean_z_lesion_temp_lvo3; mean(mean_z_lesion_temp_lvo(discriminative_features{f}, i), 1) ];
                end
                
                c = cvpartition(LOOCV_group,'resubstitution');
                [fs,history] = sequentialfs(@fs_func_svm,mean_z_lesion_temp_lvo2',LOOCV_group,'cv',c, 'direction', 'backward', 'options',opts);
                
                discriminative_features_fs = fs;
                discriminative_features_set = [ discriminative_features_set; discriminative_features_fs ];
                training_data = mean_z_lesion_temp_lvo2(discriminative_features_fs, :);
                new_data = mean_z_lesion_temp_lvo3(discriminative_features_fs); 
                
                SVMstruct = svmtrain(training_data', LOOCV_group, 'Kernel_Function', 'polynomial');
                newClass = svmclassify(SVMstruct, new_data');
                classified_histo(i) = newClass;
                SVMstruct_set(i) = { SVMstruct };
                
            end
            predict_AllMod_global_structure = (sum(strcmp(classified_histo, histo_type_fMRI)))/length(histo_type_fMRI);
            
            %% A diagnostic test
            diagnostic_test_AllMod_global_structure = [
                sum(strcmp(histo_type_fMRI(strcmp(histo_type_fMRI, 'FCDIIa'), 1), classified_histo(strcmp(histo_type_fMRI, 'FCDIIa'), 1))) sum(~(strcmp(histo_type_fMRI(strcmp(histo_type_fMRI, 'FCDIIb'), 1), classified_histo(strcmp(histo_type_fMRI, 'FCDIIb'), 1))));
                sum(~(strcmp(histo_type_fMRI(strcmp(histo_type_fMRI, 'FCDIIa'), 1), classified_histo(strcmp(histo_type_fMRI, 'FCDIIa'), 1)))) sum(strcmp(histo_type_fMRI(strcmp(histo_type_fMRI, 'FCDIIb'), 1), classified_histo(strcmp(histo_type_fMRI, 'FCDIIb'), 1)))
                ];           
            
            sum(discriminative_features_set, 1) % Vertical gradient and ALFF + ReHo
            predict_AllMod_global_case_structure = strcmp(classified_histo, histo_type_fMRI);
            
            %% ------------------------------------------------
            %%                               Histology
            %%                       FCD IIa         FCD IIb
            %% ------------------------------------------------
            %%            FCD IIa      7       |       3
            %% Prediction -------------------------------------
            %%            FCD IIb      1        |      12
            %% ------------------------------------------------
            %% sensitivity: 19/23 -> 83%
            %%        IIa : 7/8   -> 88%
            %%        IIb : 12/15 -> 80%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        end
        
    end
    
    %% 2. Use depth multi-modality features based on pre-selection
    for prediction_using_multimodal_depth_features = 1
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Feature set index
        %% mean_z_lesion
        %% 1:5   = RI of T1                      -> 2:5      p1
        %% 6:10  = PG of T1                      -> 7:9 16   p2
        %% 11:15 = TG of T1                      -> 12:15    p3
        %% 16    = PG_gw of T1                   -> 7:9 16   p2
        %% 17    = CT of T1                      -> 17       p7
        %% 18    = MC of T1                      -> 18       p8
        %% 19    = SD of T1                      -> 19       p9
        %% 20:22 = RI of T1 subcortical          -> 20:22    p4
        %% 23:25 = PG of T1 subcortical          -> 23:25    p5
        %% 26:28 = TG of T1 subcortical          -> 26:28    p6
        %% 29:33 = RI of FLAIR                   -> 30:33    p10
        %% 34:38 = PG of FLAIR                   -> 35:37 44 p11
        %% 39:43 = TG of FLAIR                   -> 40:43    p12
        %% 44    = PG_gw of FLAIR                -> 35:37 44 p11
        %% 45:47 = RI of FLAIR subcortical       -> 45:47    p13
        %% 48:50 = PG of FLAIR subcortical       -> 48:50    p14
        %% 51:53 = TG of FLAIR subcortical       -> 51:53    p15
        %% 54:58 = FA of DTI Intra               -> 56 58    p16
        %% 59:65 = FA of DT1 Subcortical         -> 60 62    p17
        %% 66:70 = MD of DTI Intra               -> 68 70    p18
        %% 71:77 = MD of DTI Subcortical         -> 72 74    p19
        %% 78    = ALFF of functional            -> 78       p20
        %% 79    = ReHo of functional            -> 79       p21
        %% 80    = DC of functional              -> 80       p22
        %% 81    = BC of functional              -> 81       p23
        %% 82 83 = Thickness of T1 perilesion    -> 82 83    p24
        %% 84 85 = Intensity of FLAIR perilesion -> 84 85    p25
        %% 86 87 = FA of DTI perilesion          -> 86 87    p26
        %% 88 89 = ALFF perilesion               -> 88 89    p27
        %% 90 91 = ReHo perilesion               -> 90 91    p28
        
        % Using t1, flair, dti and fMRI all features
        for all_modalities_together = 1
            
            classified_histo = cell(length(histo_type_fMRI), 1);
            SVMstruct_set = cell(length(histo_type_fMRI), 1);
            discriminative_features_set = [];
            opts = statset('display', 'iter');
            mean_z_lesion_final = mean_z_lesion_temp;
            keepout_vec = logical(zeros(size(mean_z_lesion_final, 1), 1));
            keepout_vec([1 6 10 11 23:28 29 34 38 39 48:53 54 55 57 59 61 63:65 66 67 69 71 73 75:77 80 81 ]) = 1;
            
            parfor i = 1 : length(histo_type_fMRI)
                
                LOOCV_sub = 1 : length(histo_type_fMRI);
                LOOCV_sub(i) = [];
                LOOCV_group = histo_type_fMRI;
                LOOCV_group(i) = [];
                
                mean_z_lesion_temp_lvo = mean_z_lesion_final;
                mean_z_lesion_temp_lvo(:, i) = [];
               
                c = cvpartition(LOOCV_group,'resubstitution');                
                [fs,history] = sequentialfs(@fs_func_svm,mean_z_lesion_temp_lvo',LOOCV_group,'cv',c, 'direction', 'forward', 'keepout', keepout_vec', 'options',opts);
                
                discriminative_features = fs;
                discriminative_features_set = [ discriminative_features_set; discriminative_features ];
                training_data = mean_z_lesion_temp_lvo(discriminative_features, :);
                new_data = mean_z_lesion_final(discriminative_features, i);        
              
                SVMstruct = svmtrain(training_data', LOOCV_group, 'Kernel_Function', 'polynomial');
                newClass = svmclassify(SVMstruct, new_data');
                classified_histo(i) = newClass;
                SVMstruct_set(i) = { SVMstruct };
                
            end
            predict_AllMod_depth = (sum(strcmp(classified_histo, histo_type_fMRI)))/length(histo_type_fMRI);            
            discriminative_features_set_AllMod_depth = discriminative_features_set;
            
            %% A diagnostic test
            diagnostic_test_AllMod_depth = [
                sum(strcmp(histo_type_fMRI(strcmp(histo_type_fMRI, 'FCDIIa'), 1), classified_histo(strcmp(histo_type_fMRI, 'FCDIIa'), 1))) sum(~(strcmp(histo_type_fMRI(strcmp(histo_type_fMRI, 'FCDIIb'), 1), classified_histo(strcmp(histo_type_fMRI, 'FCDIIb'), 1))));
                sum(~(strcmp(histo_type_fMRI(strcmp(histo_type_fMRI, 'FCDIIa'), 1), classified_histo(strcmp(histo_type_fMRI, 'FCDIIa'), 1)))) sum(strcmp(histo_type_fMRI(strcmp(histo_type_fMRI, 'FCDIIb'), 1), classified_histo(strcmp(histo_type_fMRI, 'FCDIIb'), 1)))
                ];
            
            sum(discriminative_features_set, 1) % Vertical gradient of T1 and fMRI ReHo
            predict_AllMod_depth_case = strcmp(classified_histo, histo_type_fMRI);
            
            %% ------------------------------------------------
            %%                               Histology
            %%                       FCD IIa         FCD IIb
            %% ------------------------------------------------
            %%            FCD IIa      7        |       1
            %% Prediction -------------------------------------
            %%            FCD IIb      1        |       14
            %% ------------------------------------------------
            %% sensitivity: 21/23 -> 91%
            %%        IIa : 7/8   -> 88%
            %%        IIb : 14/15 -> 93%
            
        end
        
        % Using t1, flair, dti and fMRI, and perilesion abnormalities
        for all_features_including_perilesional_ab = 1
            
            classified_histo = cell(length(histo_type_fMRI), 1);
            SVMstruct_set = cell(length(histo_type_fMRI), 1);
            discriminative_features_set = [];
            opts = statset('display', 'iter');
            mean_z_lesion_final = mean_z_lesion_temp2;
            keepout_vec = logical(zeros(size(mean_z_lesion_final, 1), 1));
            keepout_vec([1 6 10 11 23:28 29 34 38 39 48:53 54 55 57 59 61 63:65 66 67 69 71 73 75:77 80 81 ]) = 1;
                        
            parfor i = 1 : length(histo_type_fMRI)
                
                LOOCV_sub = 1 : length(histo_type_fMRI);
                LOOCV_sub(i) = [];
                LOOCV_group = histo_type_fMRI;
                LOOCV_group(i) = [];
                
                mean_z_lesion_temp_lvo = mean_z_lesion_final;
                mean_z_lesion_temp_lvo(:, i) = [];

                c = cvpartition(LOOCV_group,'resubstitution');                
                [fs,history] = sequentialfs(@fs_func_svm,mean_z_lesion_temp_lvo',LOOCV_group,'cv',c, 'direction', 'forward', 'keepout', keepout_vec', 'options',opts);

                discriminative_features = fs;
                discriminative_features_set = [ discriminative_features_set; discriminative_features ];
                training_data = mean_z_lesion_temp_lvo(discriminative_features, :);
                new_data = mean_z_lesion_final(discriminative_features, i);                
                
                SVMstruct = svmtrain(training_data', LOOCV_group, 'Kernel_Function', 'polynomial');
                newClass = svmclassify(SVMstruct, new_data');
                classified_histo(i) = newClass;
                SVMstruct_set(i) = { SVMstruct };
            end            
            predict_AllModPeri_depth = (sum(strcmp(classified_histo, histo_type_fMRI)))/length(histo_type_fMRI);
            discriminative_features_set_AllModPeri_depth = discriminative_features_set;
            
            %% A diagnostic test
            diagnostic_test_ALLModPeri_depth = [
                sum(strcmp(histo_type_fMRI(strcmp(histo_type_fMRI, 'FCDIIa'), 1), classified_histo(strcmp(histo_type_fMRI, 'FCDIIa'), 1))) sum(~(strcmp(histo_type_fMRI(strcmp(histo_type_fMRI, 'FCDIIb'), 1), classified_histo(strcmp(histo_type_fMRI, 'FCDIIb'), 1))));
                sum(~(strcmp(histo_type_fMRI(strcmp(histo_type_fMRI, 'FCDIIa'), 1), classified_histo(strcmp(histo_type_fMRI, 'FCDIIa'), 1)))) sum(strcmp(histo_type_fMRI(strcmp(histo_type_fMRI, 'FCDIIb'), 1), classified_histo(strcmp(histo_type_fMRI, 'FCDIIb'), 1)))
                ]
            
            sum(discriminative_features_set, 1) % Vertical gradient of T1, cortical intensity of FLAIR and peri-lesional cortical FA
            predict_AllModPeri_depth_case = strcmp(classified_histo, histo_type_fMRI);
            
            %% ------------------------------------------------
            %%                               Histology
            %%                       FCD IIa         FCD IIb
            %% ------------------------------------------------
            %%            FCD IIa      7        |       0
            %% Prediction -------------------------------------
            %%            FCD IIb      1        |       15
            %% ------------------------------------------------
            %% sensitivity: 22/23 -> 96%
            %%        IIa : 7/8   -> 88%
            %%        IIb : 15/15 -> 100%
            
        end
        
    end
    
    %% 3. Check Which features are most contributing for prediction
    for contributing_feature = 1
        
        feat_freq = [ sum(discriminative_features_set_AllMod_depth, 1) zeros(1, 10) ] + sum(discriminative_features_set_AllModPeri_depth, 1);
        mostcont_feat = find([ sum(discriminative_features_set_AllMod_depth, 1) zeros(1, 10) ] + sum(discriminative_features_set_AllModPeri_depth, 1));
        [a b] = sort(feat_freq(mostcont_feat), 'descend');
        mostcont_feat(b)
        
    end
    
    %% 4. Mcnemar test for comparing classifier accuracy
    for Mcnemar_test = 1
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % t1w_global vs. all_modalities_global
        mcnemar([sum(predict_T1w_global_case==1 & predict_AllMod_global_case==1) sum(predict_T1w_global_case==1 & predict_AllMod_global_case==0);
                 sum(predict_T1w_global_case==0 & predict_AllMod_global_case==1) sum(predict_T1w_global_case==0 & predict_AllMod_global_case==0) ])         
       
        % FLAIR_global_case vs. all_modalities_global
        mcnemar([sum(predict_FLAIR_global_case==1 & predict_AllMod_global_case==1) sum(predict_FLAIR_global_case==1 & predict_AllMod_global_case==0);
                 sum(predict_FLAIR_global_case==0 & predict_AllMod_global_case==1) sum(predict_FLAIR_global_case==0 & predict_AllMod_global_case==0) ])
               
        % DWI_global_case vs. all_modalities_global
        mcnemar([sum(predict_DTI_global_case==1 & predict_AllMod_global_case==1) sum(predict_DTI_global_case==1 & predict_AllMod_global_case==0);
                 sum(predict_DTI_global_case==0 & predict_AllMod_global_case==1) sum(predict_DTI_global_case==0 & predict_AllMod_global_case==0) ])

        % fMRI_global_case vs. all_modalities_global
        mcnemar([sum(predict_fMRI_global_case==1 & predict_AllMod_global_case==1) sum(predict_fMRI_global_case==1 & predict_AllMod_global_case==0);
                 sum(predict_fMRI_global_case==0 & predict_AllMod_global_case==1) sum(predict_fMRI_global_case==0 & predict_AllMod_global_case==0) ])             
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
        % t1w_global vs. all_modalities_depth
        mcnemar([sum(predict_T1w_global_case==1 & predict_AllMod_depth_case==1) sum(predict_T1w_global_case==1 & predict_AllMod_depth_case==0);
                 sum(predict_T1w_global_case==0 & predict_AllMod_depth_case==1) sum(predict_T1w_global_case==0 & predict_AllMod_depth_case==0) ])         
       
        % FLAIR_global_case vs. all_modalities_depth
        mcnemar([sum(predict_FLAIR_global_case==1 & predict_AllMod_depth_case==1) sum(predict_FLAIR_global_case==1 & predict_AllMod_depth_case==0);
                 sum(predict_FLAIR_global_case==0 & predict_AllMod_depth_case==1) sum(predict_FLAIR_global_case==0 & predict_AllMod_depth_case==0) ])
               
        % DWI_global_case vs. all_modalities_depth
        mcnemar([sum(predict_DTI_global_case==1 & predict_AllMod_depth_case==1) sum(predict_DTI_global_case==1 & predict_AllMod_depth_case==0);
                 sum(predict_DTI_global_case==0 & predict_AllMod_depth_case==1) sum(predict_DTI_global_case==0 & predict_AllMod_depth_case==0) ])   
             
        % fMRI_global_case vs. all_modalities_depth
        mcnemar([sum(predict_fMRI_global_case==1 & predict_AllMod_depth_case==1) sum(predict_fMRI_global_case==1 & predict_AllMod_depth_case==0);
                 sum(predict_fMRI_global_case==0 & predict_AllMod_depth_case==1) sum(predict_fMRI_global_case==0 & predict_AllMod_depth_case==0) ]) 
       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % all_modalities_global vs. all_modalities_depth
        mcnemar([sum(predict_AllMod_global_case==1 & predict_AllMod_depth_case==1) sum(predict_AllMod_global_case==1 & predict_AllMod_depth_case==0);
                 sum(predict_AllMod_global_case==0 & predict_AllMod_depth_case==1) sum(predict_AllMod_global_case==0 & predict_AllMod_depth_case==0) ])
             
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
        % t1w_global vs. all_modalities_depth + peri
        mcnemar([sum(predict_T1w_global_case==1 & predict_AllModPeri_depth_case==1) sum(predict_T1w_global_case==1 & predict_AllModPeri_depth_case==0);
                 sum(predict_T1w_global_case==0 & predict_AllModPeri_depth_case==1) sum(predict_T1w_global_case==0 & predict_AllModPeri_depth_case==0) ])         
       
        % FLAIR_global_case vs. all_modalities_depth  + peri
        mcnemar([sum(predict_FLAIR_global_case==1 & predict_AllModPeri_depth_case==1) sum(predict_FLAIR_global_case==1 & predict_AllMod_depth_case==0);
                 sum(predict_FLAIR_global_case==0 & predict_AllModPeri_depth_case==1) sum(predict_FLAIR_global_case==0 & predict_AllModPeri_depth_case==0) ])
               
        % DWI_global_case vs. all_modalities_depth  + peri
        mcnemar([sum(predict_DTI_global_case==1 & predict_AllModPeri_depth_case==1) sum(predict_DTI_global_case==1 & predict_AllModPeri_depth_case==0);
                 sum(predict_DTI_global_case==0 & predict_AllModPeri_depth_case==1) sum(predict_DTI_global_case==0 & predict_AllModPeri_depth_case==0) ])   
             
        % fMRI_global_case vs. all_modalities_depth + peri
        mcnemar([sum(predict_fMRI_global_case==1 & predict_AllModPeri_depth_case==1) sum(predict_fMRI_global_case==1 & predict_AllModPeri_depth_case==0);
                 sum(predict_fMRI_global_case==0 & predict_AllModPeri_depth_case==1) sum(predict_fMRI_global_case==0 & predict_AllModPeri_depth_case==0) ]) 
       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % all_modalities_global vs. all_modalities_depth + peri
        mcnemar([sum(predict_AllMod_global_case==1 & predict_AllModPeri_depth_case==1) sum(predict_AllMod_global_case==1 & predict_AllModPeri_depth_case==0);
                 sum(predict_AllMod_global_case==0 & predict_AllModPeri_depth_case==1) sum(predict_AllMod_global_case==0 & predict_AllModPeri_depth_case==0) ])    
             
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % all_modalities_depth vs. all_modalities_depth + peri
        mcnemar([sum(predict_AllMod_depth_case==1 & predict_AllModPeri_depth_case==1) sum(predict_AllMod_depth_case==1 & predict_AllModPeri_depth_case==0);
                 sum(predict_AllMod_depth_case==0 & predict_AllModPeri_depth_case==1) sum(predict_AllMod_depth_case==0 & predict_AllModPeri_depth_case==0) ])
        
    end
    
    %% 5. Do post-hoc analysis for a missing case
    for missing_case = 1
        
        % This case, most features seem to fit well the pattern of Type IIA
        % but cortical and subcortical FA and T1 vertical gradient at the
        % GM-WM 
        idx = 1:91;
        idx = idx([ 2:5 7:9 16 12:15 17 18 19 20:22 30:33 35:37 44 40:43 45:47 56 58 60 62 68 70 72 74 78 79 82:91 ]);
%         featvec_missing = mean_z_lesion_temp2([ 2:5 7:9 16 12:15 17 18 19 20:22 30:33 35:37 44 40:43 45:47 56 58 60 62 68 70 72 74 78 79 82:91 ], ~predict_AllModPeri_depth_case);
        featvec_missing = mean_z_lesion_temp2([ 2:5 7:9 16 12:15 17 18 19 20:22 30:33 35:37 44 40:43 45:47 56 58 60 62 68 70 72 74 78 79 82:91 ], ~predict_AllMod_depth_case);
        [ idx' featvec_missing ]
        histo_type_fMRI(~predict_AllModPeri_depth_case)
        featvec_typeIIA = mean(mean_z_lesion_temp2([ 2:5 7:9 16 12:15 17 18 19 20:22 30:33 35:37 44 40:43 45:47 56 58 60 62 68 70 72 74 78 79 82:91 ], strcmp(histo_type_fMRI, 'FCDIIa')), 2)
        featvec_typeIIB = mean(mean_z_lesion_temp2([ 2:5 7:9 16 12:15 17 18 19 20:22 30:33 35:37 44 40:43 45:47 56 58 60 62 68 70 72 74 78 79 82:91 ], strcmp(histo_type_fMRI, 'FCDIIb')), 2)        
        
        %% missing two cases -> 071
        abs(featvec_missing(:, 1) - featvec_typeIIA)/sum(abs(featvec_missing(:, 1) - featvec_typeIIA)) % which feature contributed most to be different from type IIA?
        abs(featvec_missing(:, 1) - featvec_typeIIB)/sum(abs(featvec_missing(:, 1) - featvec_typeIIB)) % which feature contributed most to be different from type IIB?
        sqrt(sum((featvec_missing(:, 1) - featvec_typeIIA).^2)) % How much different it is from type IIA?
        sqrt(sum((featvec_missing(:, 1) - featvec_typeIIB).^2)) % How much different it is from type IIB?
        
    end    
    
    %% 6. Visualize the separation
    for visualizatoin_using_MDS = 1
        
        discriminative_features = find(sum(discriminative_features_set_AllModPeri_depth, 1) > 0);                
        % ================================
        % ===== Showing the results ======
        % ================================
        case_num_pat_temp = case_num_pat_fMRI; % case_num_pat_temp(17) =[];
        
        % Assign color for each class
        colorList = [ 1 0 0; 0 0 0 ];
        
        % multidimensional scaling
        for mds = 1
            
            % true (ground truth) class
            N = length(histo_type_fMRI);
            trueClassIndex = zeros(N,1);
            trueClassIndex(strcmp(histo_type_fMRI, 'FCDIIa')) = 1;
            trueClassIndex(strcmp(histo_type_fMRI, 'FCDIIb')) = 2; 
            colorTrueClass = colorList(trueClassIndex,:);
            % result Class
            resultClassIndex = zeros(N,1);
            resultClassIndex(strcmp(classified_histo, 'FCDIIa')) = 1;
            resultClassIndex(strcmp(classified_histo, 'FCDIIb')) = 2; 
            colorResultClass = colorList(resultClassIndex,:);
            
            % Reduce the dimension from 15D to 2D
            data = mean_z_lesion_temp2(discriminative_features, :); 
            distanceMatrix = pdist(data','euclidean');
            [ newCoor stress disparities ] = mdscale(distanceMatrix,2);
            
            x = newCoor(:,1); x(11) = 0.3;
            y = newCoor(:,2); y(12) = 1.8; y(21) = -1.9;
            patchSize = 30; 
            colorTrueClassPlot = colorTrueClass;
            colorResultClassPlot = colorResultClass;
            figure;
            position = get(gcf, 'Position'); position = [ 2200 600 1000 600 ]; set(gcf, 'Position', position);
            S1 = subplot(1, 2, 1); scatter(x,y,patchSize,colorTrueClassPlot,'filled');  xlim([-3.5 3.5]); ylim([-3.5 3.5]); title('actual grouping');
            S2 = subplot(1, 2, 2); scatter(x,y,patchSize,colorResultClassPlot,'filled');  xlim([-3.5 3.5]); ylim([-3.5 3.5]); title('classificaition');
            set( get(S1, 'XLabel'), 'String', 'x component'); set( get(S2, 'XLabel'), 'String', 'x componenet');
            set( get(S1, 'YLabel'), 'String', 'y component');
            [ax, h] = suplabel('Multidimensional scaling', 't'); set(h,'FontSize',15);
            
        end
        
        for mds_final_1 = 1
            
            % true (ground truth) class
            N = length(histo_type_fMRI);
            trueClassIndex = zeros(N,1);
            trueClassIndex(strcmp(histo_type_fMRI, 'FCDIIa')) = 1;
            trueClassIndex(strcmp(histo_type_fMRI, 'FCDIIb')) = 2; 
            colorTrueClass = colorList(trueClassIndex,:);
            % result Class
            resultClassIndex = zeros(N,1);
            resultClassIndex(strcmp(classified_histo, 'FCDIIa')) = 1;
            resultClassIndex(strcmp(classified_histo, 'FCDIIb')) = 2; 
            colorResultClass = colorList(resultClassIndex,:);
            
            misclassified = find(resultClassIndex-trueClassIndex);
            
            % Reduce the dimension from 15D to 2D
            data = mean_z_lesion_temp2(discriminative_features, :);             
            distanceMatrix = pdist(data','euclidean');
            [ newCoor stress disparities ] = mdscale(distanceMatrix,2);
            
            x = newCoor(:,1); x(11) = 0.3;
            y = newCoor(:,2); y(12) = 1.8; y(21) = -1.9;
            patchSize = 30; 
            colorTrueClassPlot = colorTrueClass;
            colorResultClassPlot = colorResultClass;
            figure; S1 = axis; hold on;      
            scatter(x,y,patchSize,colorTrueClassPlot,'filled'); 
            scatter(x,y,patchSize*2,colorResultClass); 
            scatter(x(misclassified),y(misclassified),patchSize*2,colorResultClass(misclassified), 'x');
            xlim([-3 3.5]); ylim([-2.5 2]); title('Multidimensional scaling');
            xlabel('x component'); ylabel('y component');  
            
            % legend
            scatter(-6,7,patchSize,colorList(1, :),'filled'); text(-5.5, 7, 'IIA');
            scatter(-6,6.3,  patchSize,colorList(2, :),'filled'); text(-5.5, 6.3, 'IIB');
            scatter(-6,5.4,patchSize*1.7,colorList(2, :)); text(-5.5, 5.4, 'correctly classified');
            scatter(-6,4.5,  patchSize*1.7,colorList(2, :)); scatter(-6,4.5,  patchSize*1.7,colorList(2, :), 'x'); text(-5.5, 4.5, 'misclassifed');
            
            scatter(x,y,patchSize*2,colorResultClass); 
            scatter(x(misclassified),y(misclassified),patchSize*2,colorResultClass(misclassified), 'x');
            
            if(save_fig)
                export_fig([OUTPATH '/prediction_mds_summary_final_new' ], '-m4', '-png'); close(gcf);
            end
            
        end
        
        discriminative_features = find(sum(discriminative_features_set_AllMod_depth, 1) > 0);        
        % ================================
        % ===== Showing the results ======
        % ================================
        case_num_pat_temp = case_num_pat_fMRI; % case_num_pat_temp(17) =[];
        
        % Assign color for each class
        colorList = [ 1 0 0; 0 0 0 ];
        
        for mds_final_2 = 1
            
            % true (ground truth) class
            N = length(histo_type_fMRI);
            trueClassIndex = zeros(N,1);
            trueClassIndex(strcmp(histo_type_fMRI, 'FCDIIa')) = 1;
            trueClassIndex(strcmp(histo_type_fMRI, 'FCDIIb')) = 2; 
            colorTrueClass = colorList(trueClassIndex,:);
            % result Class
            resultClassIndex = zeros(N,1);
            resultClassIndex(strcmp(classified_histo, 'FCDIIa')) = 1;
            resultClassIndex(strcmp(classified_histo, 'FCDIIb')) = 2; 
            colorResultClass = colorList(resultClassIndex,:);
            
            misclassified = find(resultClassIndex-trueClassIndex);
            
            % Reduce the dimension from 15D to 2D
            data = mean_z_lesion_temp(discriminative_features, :); 
            distanceMatrix = pdist(data','euclidean');
            [ newCoor stress disparities ] = mdscale(distanceMatrix,2);
            
            x = newCoor(:,1); x(15) = -0.9; x(22) = -0.3;
            y = newCoor(:,2); y(3) = -1.5;
            patchSize = 30; 
            colorTrueClassPlot = colorTrueClass;
            colorResultClassPlot = colorResultClass;
            figure; S1 = axis; hold on;      
            scatter(x,y,patchSize,colorTrueClassPlot,'filled'); 
            scatter(x,y,patchSize*2,colorResultClass); 
            scatter(x(misclassified),y(misclassified),patchSize*2,colorResultClass(misclassified), 'x');
            xlim([-3 3]); ylim([-2.5 2.5]); title('Multidimensional scaling');
            xlabel('x component'); ylabel('y component');  
            
            % legend
            scatter(-6,7,patchSize,colorList(1, :),'filled'); text(-5.5, 7, 'IIA');
            scatter(-6,6.3,  patchSize,colorList(2, :),'filled'); text(-5.5, 6.3, 'IIB');
            scatter(-6,5.4,patchSize*1.7,colorList(2, :)); text(-5.5, 5.4, 'correctly classified');
            scatter(-6,4.5,  patchSize*1.7,colorList(2, :)); scatter(-6,4.5,  patchSize*1.7,colorList(2, :), 'x'); text(-5.5, 4.5, 'misclassifed');
            
            scatter(x,y,patchSize*2,colorResultClass); 
            scatter(x(misclassified),y(misclassified),patchSize*2,colorResultClass(misclassified), 'x');
            
            if(save_fig)
                export_fig([OUTPATH '/prediction_mds_summary_final_new2' ], '-m4', '-png'); close(gcf);
            end
            
        end        
        
    end
    
    %% 7. Permutation test for checking whether it is over the chance level
    for statistical_test_using_all_features_with_perilesional_ab = 1
                
        permutation = 1000;
        prediction_accuracy_rand = zeros(permutation, 1);
        discriminative_features_set = [];        
        mean_z_lesion_final = mean_z_lesion_temp2;
        keepout_vec = logical(zeros(size(mean_z_lesion_final, 1), 1));
        keepout_vec([1 6 10 11 23:28 29 34 38 39 48:53 54 55 57 59 61 63:65 66 67 69 71 73 75:77 80 81 ]) = 1;
        
        for p = 1 : permutation
            
            rand_idx = randperm(length(histo_type_fMRI));
            histo_type_rand = histo_type_fMRI(rand_idx);   
            
            classified_histo = cell(length(histo_type_rand), 1);
            
            parfor i = 1 : length(histo_type_rand)
                
                LOOCV_sub = 1 : length(histo_type_rand);
                LOOCV_sub(i) = [];
                LOOCV_group = histo_type_rand;
                LOOCV_group(i) = [];
                
                mean_z_lesion_temp_lvo = mean_z_lesion_final;
                mean_z_lesion_temp_lvo(:, i) = [];

                c = cvpartition(LOOCV_group,'resubstitution');                
                [fs,history] = sequentialfs(@fs_func_svm,mean_z_lesion_temp_lvo',LOOCV_group,'cv',c, 'direction', 'forward', 'keepout', keepout_vec');

                discriminative_features = fs;
                discriminative_features_set = [ discriminative_features_set; discriminative_features ];
                training_data = mean_z_lesion_temp_lvo(discriminative_features, :);
                new_data = mean_z_lesion_final(discriminative_features, i);                
                
                SVMstruct = svmtrain(training_data', LOOCV_group, 'Kernel_Function', 'polynomial');
                newClass = svmclassify(SVMstruct, new_data');
                classified_histo(i) = newClass;

            end
            
            disp([ 'iter ' num2str(p) ': ' num2str((sum(strcmp(classified_histo, histo_type_rand)))/length(histo_type_rand)) ]);
            prediction_accuracy_rand(p) = (sum(strcmp(classified_histo, histo_type_rand)))/length(histo_type_rand);
            
        end
        
    end
    
    (sum(predict_AllModPeri_depth < sort(prediction_accuracy_rand)) + 1)/permutation
    mean(prediction_accuracy_rand)
    std(prediction_accuracy_rand)

end