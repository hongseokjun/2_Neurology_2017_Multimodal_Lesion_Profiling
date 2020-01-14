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

    load('zscore_database_T1_FLAIR_cont.mat');
    mean_z_lesion_T1_FLAIR_cont = mean_z_lesion_cont;
    load('zscore_database_DTI_cont.mat');
    mean_z_lesion_DTI_cont = mean_z_lesion_cont;
    load('zscore_database_fMRI_cont.mat');
    mean_z_lesion_fMRI_cont = mean_z_lesion_cont;
    
    mean_z_lesion_T1_FLAIR_cont_org = mean_z_lesion_T1_FLAIR_cont;
    mean_z_lesion_DTI_cont_org = mean_z_lesion_DTI_cont;
    mean_z_lesion_fMRI_cont_org = mean_z_lesion_fMRI_cont;
    
    %% perilesional profiling
    Kernel = 2;    
    load(['pl_feature_set_sm_' num2str(Kernel) '_T1_FLAIR.mat']);
    pl_feature_set_T1_FLAIR = pl_feature_set;
    load(['pl_feature_set_sm_' num2str(Kernel) '_DTI.mat']);
    pl_feature_set_DTI = pl_feature_set;
    load(['pl_feature_set_sm_' num2str(Kernel) '_fMRI.mat']);
    pl_feature_set_fMRI = pl_feature_set;

    load(['pl_feature_set_sm_' num2str(Kernel) '_T1_FLAIR_cont.mat']);
    pl_feature_set_T1_FLAIR_cont = pl_feature_set;
    load(['pl_feature_set_sm_' num2str(Kernel) '_DTI_cont.mat']);
    pl_feature_set_DTI_cont = pl_feature_set;
    load(['pl_feature_set_sm_' num2str(Kernel) '_fMRI_cont.mat']);
    pl_feature_set_fMRI_cont = pl_feature_set;
    
    pl_feature_set_T1_FLAIR_cont_org = pl_feature_set_T1_FLAIR_cont;
    pl_feature_set_DTI_cont_org = pl_feature_set_DTI_cont;
    pl_feature_set_fMRI_cont_org = pl_feature_set_fMRI_cont;

    %% statistical test and permutation result
    load('zscore_database_perilesional_T1_FLAIR_DTI_fMRI_prediction.mat');
    load('permutation_test_classifier.mat');
    
end

%% Read demograpic data and set up some parameters and variables
for read_demodata_setup_params = 1
    
    BASEDIR = '
    OUTPATH = '/local_raid/seokjun/01_project/04_IntracorticalAnalysis/03_result/prediction_histology_by_MRI/';
    Cases_cont                = '/local_raid/seokjun/01_project/04_IntracorticalAnalysis/01_analysis/Demographic_data_control.txt';
    Cases_cont_dti            = '/local_raid/seokjun/01_project/04_IntracorticalAnalysis/01_analysis/Demographic_data_control_DTI.txt';
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

    %% read controls data (t1, flair)
    fid = fopen(Cases_cont);
    demo = textscan(fid, '%s%f%s', 'Delimiter', ',', 'CollectOutput', 1);
    case_num_cont = demo{1};    
    age_cont = demo{2};
    gender_cont = demo{3};
    case_num_cont_org = case_num_cont;
    fclose(fid);

    %% read controls data (dti)
    fid = fopen(Cases_cont_dti);
    demo = textscan(fid, '%s%f%s', 'Delimiter', ',', 'CollectOutput', 1);
    case_num_cont_dti = demo{1};
    case_num_cont_dti_org = case_num_cont_dti;
    fclose(fid);
    
    %% read controls data (fMRI)
    fid = fopen(Cases_cont_fMRI);
    demo = textscan(fid, '%s%f%s%s%d', 'Delimiter', ',', 'CollectOutput', 1);
    case_num_cont_fMRI = demo{1};    
    case_num_cont_fMRI(strcmp(case_num_cont_fMRI, '322_2')) = {'322_1'};    
    case_num_cont_fMRI_org = case_num_cont_fMRI;
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
    
    %% 1. Make the feature vectors, matching the T1, FLAIR and DTI zscore
    for featureset_build = 1
        
        for control_setup = 1
            
            case_num_cont = case_num_cont_org;
            case_num_cont_dti = case_num_cont_dti_org;
            case_num_cont_fMRI = case_num_cont_fMRI_org;
            pl_feature_set_T1_FLAIR_cont = pl_feature_set_T1_FLAIR_cont_org;
            pl_feature_set_DTI_cont = pl_feature_set_DTI_cont_org;
            pl_feature_set_fMRI_cont = pl_feature_set_fMRI_cont_org;
            mean_z_lesion_T1_FLAIR_cont = mean_z_lesion_T1_FLAIR_cont_org;
            mean_z_lesion_DTI_cont = mean_z_lesion_DTI_cont_org;
            mean_z_lesion_fMRI_cont = mean_z_lesion_fMRI_cont_org;
            
            %% fMRI
            excluded_cases = { '306_1', '313_1', '322_1', '342_1', '343_1' };
            [C, ia1, ib1] = intersect(case_num_cont_fMRI, excluded_cases);
            [C, iap, ibp] = intersect(case_num_pat, case_num_pat_fMRI);
            case_num_cont_fMRI(ia1) = [];
            mean_z_lesion_fMRI_cont(:, :, ia1) = [];
            for i = 1 : size(pl_feature_set_fMRI_cont, 2)
                pl_feature_set_fMRI_cont(i).data = pl_feature_set_fMRI_cont(i).data(:, :, :, iap);
            end
            
            %% T1/FLAIR
            excluded_cases = { '306_1', '313_1', '322_1' };
            [C, ia1, ib1] = intersect(case_num_cont, excluded_cases);
            case_num_cont(ia1) = [];
            [C, ia1, ib1] = intersect(case_num_cont, case_num_cont_fMRI);
            case_num_cont = case_num_cont(ia1);
            mean_z_lesion_T1_FLAIR_cont = mean_z_lesion_T1_FLAIR_cont(:, :, ia1);
            for i = 1 : size(pl_feature_set_T1_FLAIR_cont, 2)
                pl_feature_set_T1_FLAIR_cont(i).data = pl_feature_set_T1_FLAIR_cont(i).data(:, ia1, :, :);
            end
            
            %% DTI
            excluded_cases = { '306_1', '313_1', '322_1' };
            [C, ia1, ib1] = intersect(case_num_cont_dti, excluded_cases);
            case_num_cont_dti(ia1) = [];
            [C, ia1, ib1] = intersect(case_num_cont_dti, case_num_cont_fMRI);
            case_num_cont_dti = case_num_cont_dti(ia1);
            mean_z_lesion_DTI_cont = mean_z_lesion_DTI_cont(:, :, ia1);
            for i = 1 : size(pl_feature_set_DTI_cont, 2)
                pl_feature_set_DTI_cont(i).data = pl_feature_set_DTI_cont(i).data(:, ia1, :, :);
            end
            
            mean_z_perilesion_T1_cont_final    = mean(pl_feature_set_T1_FLAIR_cont(6).data, 4)'; mean_z_perilesion_T1_cont_final = mean(mean_z_perilesion_T1_cont_final(:, 1:3), 2)';
            mean_z_perilesion_FLAIR_cont_final = mean(squeeze(mean(pl_feature_set_T1_FLAIR_cont(9).data(:, :, 2:5, :), 3)), 3)'; mean_z_perilesion_FLAIR_cont_final = mean(mean_z_perilesion_FLAIR_cont_final(:, 1:6), 2)';
            mean_z_perilesion_DTI_cont_final   = mean(squeeze(mean(pl_feature_set_DTI_cont(1).data(:, :, 2:5, :), 3)), 3)'; mean_z_perilesion_DTI_cont_final = mean(mean_z_perilesion_DTI_cont_final(:, 1:5), 2)';
            mean_z_perilesion_fMRI_cont_final  = mean(pl_feature_set_fMRI_cont(2).data, 4)'; mean_z_perilesion_fMRI_cont_final = mean(mean_z_perilesion_fMRI_cont_final(:, 1:14), 2)';
            
            mean_z_lesion_cont      = [ squeeze(mean(mean_z_lesion_T1_FLAIR_cont, 2)); squeeze(mean(mean_z_lesion_DTI_cont, 2)); squeeze(mean(mean_z_lesion_fMRI_cont, 2)); mean_z_perilesion_T1_cont_final; mean_z_perilesion_FLAIR_cont_final; mean_z_perilesion_DTI_cont_final; mean_z_perilesion_fMRI_cont_final ];
            mean_z_lesion_cont_temp = [ squeeze(mean(mean_z_lesion_T1_FLAIR_cont, 2)); squeeze(mean(mean_z_lesion_DTI_cont, 2)); squeeze(mean(mean_z_lesion_fMRI_cont, 2)); ];
            
            peri_lesional_dist      =  [ 2:2:4 ];
            peri_lesional_dist_fMRI =  [ 2:4:8 ];
            mean_z_lesion_cont_temp2 = [ squeeze(mean(mean_z_lesion_T1_FLAIR_cont, 2)); squeeze(mean(mean_z_lesion_DTI_cont, 2)); squeeze(mean(mean_z_lesion_fMRI_cont, 2));
                                         mean(pl_feature_set_T1_FLAIR_cont(6).data(peri_lesional_dist, :, :, :), 4);                                                                        %% T1 cortical thickness
                                         mean(squeeze(mean(pl_feature_set_T1_FLAIR_cont(9).data(peri_lesional_dist, :, 2:5, :), 3)), 3)                                                     %% FLAIR cortical intentisy
                                         mean(squeeze(mean(pl_feature_set_DTI_cont(1).data(peri_lesional_dist, :, 2:5, :), 3)), 3);                                                         %% DTI cortical FA
                                         mean(pl_feature_set_fMRI_cont(1).data(peri_lesional_dist_fMRI, :, :, :), 4);                                                                       %% fMRI ALFF
                                         mean(pl_feature_set_fMRI_cont(2).data(peri_lesional_dist_fMRI, :, :, :), 4);                                                                       %% fMRI ReHo
                                       ];
            
        end

        for patient_setup = 1
            
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
            
            mean_z_lesion_FCD      = [ mean_z_lesion_T1_FLAIR(:, ia1); mean_z_lesion_DTI(:, ia2); mean_z_lesion_fMRI; mean(mean_z_perilesion_T1_final, 1); mean(mean_z_perilesion_FLAIR_final, 1); mean(mean_z_perilesion_DTI_final, 1); mean(mean_z_perilesion_fMRI_final, 1) ];
            mean_z_lesion_FCD_temp = [ mean_z_lesion_T1_FLAIR(:, ia1); mean_z_lesion_DTI(:, ia2); mean_z_lesion_fMRI; ];
            
            peri_lesional_dist      =  [ 2:2:4 ];
            peri_lesional_dist_fMRI =  [ 2:4:8 ];
            mean_z_lesion_FCD_temp2 = [ mean_z_lesion_T1_FLAIR(:, ia1); mean_z_lesion_DTI(:, ia2); mean_z_lesion_fMRI;
                                    pl_feature_set_T1_FLAIR(6).data(peri_lesional_dist, ia1);                                                                                          %% T1 cortical thickness
                                    mean(pl_feature_set_T1_FLAIR(9).data(peri_lesional_dist, ia1, 2:5), 3);                                                                            %% FLAIR cortical intentisy
                                    mean(pl_feature_set_DTI(1).data(peri_lesional_dist, ia2, [3 5]), 3);                                                                               %% DTI cortical FA
                                    pl_feature_set_fMRI(1).data(peri_lesional_dist_fMRI, :);                                                                                           %% fMRI ALFF
                                    pl_feature_set_fMRI(2).data(peri_lesional_dist_fMRI, :);                                                                                           %% fMRI ReHo
                                  ];
            
        end
        
        mean_z_lesion        = [ mean_z_lesion_cont mean_z_lesion_FCD ];
        mean_z_lesion_temp   = [ mean_z_lesion_cont_temp mean_z_lesion_FCD_temp ];
        mean_z_lesion_temp2  = [ mean_z_lesion_cont_temp2 mean_z_lesion_FCD_temp2 ];
        
        group_class = cell(size(mean_z_lesion_cont, 2)+size(mean_z_lesion_FCD, 2), 1);
        group_class(1:size(mean_z_lesion_cont, 2)) = { 'control' };
        group_class(size(mean_z_lesion_cont, 2)+1:end) = histo_type_fMRI;
        
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
        
        % Using t1, flair, dti and fMRI, and perilesion abnormalities
        for all_features_including_perilesional_ab = 1
            
            delete(gcp);
            
            parpool(24);
            
            % one vs. one 
            opts = [];
            mean_z_lesion_final = mean_z_lesion_temp2; % mean_z_lesion_final(79, 6) = 0.7;
            keepout_vec = logical(zeros(size(mean_z_lesion_final, 1), 1));
            keepout_vec([ 1 6 10 11 23:28 29 34 38 39 48:53 54 55 57 59 61 63:65 66 67 69 71 73 75:77 80 81 ]) = 1;
            
            predict_best = 0;
            discriminative_features_set_best = [];
            SVMstruct_set_best = [];
            
            permutation = 100;
            predict_AllModPeri_depth = zeros(permutation, 1);
            discriminative_features_set_AllModPeri_depth = zeros(length(group_class), size(mean_z_lesion_final, 1), permutation);
            for p = 1 : permutation
                
                classified_group = group_class;                                          
                SVMstruct_set = cell(length(classified_group), 1);
                discriminative_features_set = cell(length(group_class), 1);
                
                tic
                parfor i = 1 : length(group_class)
                    
                    LOOCV_sub = 1 : length(group_class);
                    LOOCV_sub(i) = [];
                    LOOCV_group = group_class;
                    LOOCV_group(i) = [];
                    
                    mean_z_lesion_temp_lvo = mean_z_lesion_final;
                    mean_z_lesion_temp_lvo(:, i) = [];
                    
                    c = cvpartition(LOOCV_group,'resubstitution');
                    [fs,history] = sequentialfs(@fs_func_multiclass_svm,mean_z_lesion_temp_lvo',LOOCV_group,'cv',c, 'direction', 'forward', 'keepout', keepout_vec', 'options',opts);
                    
                    discriminative_features = fs;
                    discriminative_features_set(i) = { discriminative_features };
%                     discriminative_features = logical(discriminative_features_set_AllModPeri_depth(i, :, p));
                    
                    training_data = mean_z_lesion_temp_lvo(discriminative_features, :);
                    new_data = mean_z_lesion_final(discriminative_features, i);
                    
                    t = templateSVM('KernelFunction', 'polynomial');
                    classOrder = unique(group_class);
                    
                    Mdl = fitcecoc(training_data',LOOCV_group, 'Learners', t, 'ClassNames', classOrder, 'Verbose',2);
                    newClass = predict(Mdl, new_data');
                    classified_group(i) = newClass;
                    SVMstruct_set(i) = { Mdl };
                    
                end
                toc
                
                predict_AllModPeri_depth(p) = (sum(strcmp(classified_group, group_class)))/length(group_class)
                discriminative_features_set_AllModPeri_depth(:, :, p) = cell2mat(discriminative_features_set);
                
            end           
            
            delete(gcp);
            
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
        
        p = 56;
        discriminative_features = find(mean(discriminative_features_set_AllModPeri_depth(:, :, p), 1) > 0.1);
        % ================================
        % ===== Showing the results ======
        % ================================
        case_num_pat_temp = case_num_pat_fMRI; % case_num_pat_temp(17) =[];
        
        % Assign color for each class
        colorList = [ 0 0.8 0; 0.5 0.5 0.5; 0 0 0 ];
        
        % multidimensional scaling
        for mds_1 = 1
            
            % true (ground truth) class
            N = length(group_class);
            trueClassIndex = zeros(N,1);
            trueClassIndex(strcmp(group_class, 'control')) = 1;
            trueClassIndex(strcmp(group_class, 'FCDIIa')) = 2;
            trueClassIndex(strcmp(group_class, 'FCDIIb')) = 3;
            
            colorTrueClass = colorList(trueClassIndex,:);
            % result Class
            resultClassIndex = zeros(N,1);
            resultClassIndex(strcmp(classified_group, 'control')) = 1;
            resultClassIndex(strcmp(classified_group, 'FCDIIa')) = 2;
            resultClassIndex(strcmp(classified_group, 'FCDIIb')) = 3; 
            colorResultClass = colorList(resultClassIndex,:);
            
            % Reduce the dimension from 15D to 2D
            data = mean_z_lesion_temp2(discriminative_features, :); 
            distanceMatrix = pdist(data','euclidean');
            [ newCoor stress disparities ] = mdscale(distanceMatrix,2);
            
            x = newCoor(:,1); 
            y = newCoor(:,2); 
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
                       
        for mds_2 = 1
            
            N = length(group_class);
            trueClassIndex = zeros(N,1);
            trueClassIndex(strcmp(group_class, 'control')) = 1;
            trueClassIndex(strcmp(group_class, 'FCDIIa')) = 2;
            trueClassIndex(strcmp(group_class, 'FCDIIb')) = 3;
            
            colorTrueClass = colorList(trueClassIndex,:);
            % result Class
            resultClassIndex = zeros(N,1);
            resultClassIndex(strcmp(classified_group, 'control')) = 1;
            resultClassIndex(strcmp(classified_group, 'FCDIIa')) = 2;
            resultClassIndex(strcmp(classified_group, 'FCDIIb')) = 3; 
            colorResultClass = colorList(resultClassIndex,:);
            
            % Reduce the dimension from 15D to 2D
            data = mean_z_lesion_temp2(discriminative_features, :); 
            distanceMatrix = pdist(data','euclidean');
            [ newCoor stress disparities ] = mdscale(distanceMatrix,2);
            
            x = newCoor(:,1); 
            y = newCoor(:,2); 
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
            scatter(-6,7,patchSize,colorList(1, :),'filled'); text(-5.5, 7.7, 'Controls');
            scatter(-6,7,patchSize,colorList(2, :),'filled'); text(-5.5, 7, 'IIA');
            scatter(-6,6.3,  patchSize,colorList(3, :),'filled'); text(-5.5, 6.3, 'IIB');
            scatter(-6,5.4,patchSize*1.7,colorList(2, :)); text(-5.5, 5.4, 'correctly classified');
            scatter(-6,4.5,  patchSize*1.7,colorList(2, :)); scatter(-6,4.5,  patchSize*1.7,colorList(2, :), 'x'); text(-5.5, 4.5, 'misclassifed');
            
            scatter(x,y,patchSize*2,colorResultClass); 
            scatter(x(misclassified),y(misclassified),patchSize*2,colorResultClass(misclassified), 'x');
            
            if(save_fig)
                export_fig([OUTPATH '/prediction_multiclass_mds_summary_final' ], '-m4', '-png'); close(gcf);
            end
            
        end        
        
    end
    
    %% 7. Permutation test for checking whether it is over the chance level
    for statistical_test_using_all_features_with_perilesional_ab = 1
        
        delete(gcp);
        parpool(24);
                
        permutation = 1000;
        prediction_accuracy_rand = zeros(permutation, 1);        
        opts = [];
        mean_z_lesion_final = mean_z_lesion_temp2;
        keepout_vec = logical(zeros(size(mean_z_lesion_final, 1), 1));
        keepout_vec([ 1 6 10 11 23:28 29 34 38 39 48:53 54 55 57 59 61 63:65 66 67 69 71 73 75:77 80 81 ]) = 1;
        
        % p = 0.001; 23% - 62%
        for p = 1 : permutation
            
            rand_idx = randperm(length(group_class));
            histo_type_rand = group_class(rand_idx);   
            
            classified_group = cell(length(histo_type_rand), 1);
            
            parfor i = 1 : length(histo_type_rand)
                
                LOOCV_sub = 1 : length(histo_type_rand);
                LOOCV_sub(i) = [];
                LOOCV_group = histo_type_rand;
                LOOCV_group(i) = [];
                
                mean_z_lesion_temp_lvo = mean_z_lesion_final;
                mean_z_lesion_temp_lvo(:, i) = [];

                c = cvpartition(LOOCV_group,'resubstitution');                
                [fs,history] = sequentialfs(@fs_func_multiclass_svm,mean_z_lesion_temp_lvo',LOOCV_group,'cv',c, 'direction', 'forward', 'keepout', keepout_vec');

                discriminative_features = fs;                
                training_data = mean_z_lesion_temp_lvo(discriminative_features, :);
                new_data = mean_z_lesion_final(discriminative_features, i);                
               
                Mdl = fitcecoc(training_data',LOOCV_group, 'Learners', t, 'ClassNames', classOrder, 'Verbose',2);
                newClass = predict(Mdl, new_data');
                classified_group(i) = newClass;
                SVMstruct_set(i) = { Mdl };

            end
            
            disp([ 'iter ' num2str(p) ': ' num2str((sum(strcmp(classified_group, histo_type_rand)))/length(histo_type_rand)) ]);
            prediction_accuracy_rand(p) = (sum(strcmp(classified_group, histo_type_rand)))/length(histo_type_rand);
            
        end
        
        delete(gcp);
        
    end
    
    (sum(max(predict_AllModPeri_depth) < sort(prediction_accuracy_rand)) + 1)/permutation
    mean(prediction_accuracy_rand)
    std(prediction_accuracy_rand)

end