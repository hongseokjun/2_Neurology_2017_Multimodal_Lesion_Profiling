clear all
close all
addpath('/data/noel/noel2/local/matlab/surfstat_chicago');
addpath('/local_raid/seokjun/03_downloads/LSSVMlabv1_8_R2009b_R2011a');

%% Read pre-calculated zscore of T1 and FLAIR in the lesion
for read_matfiles = 1
    
    kfold = 5; % 10;
    
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
    if kfold == 10
        load('Prediction_Histology_Result_kfold_multiclass_trial3_same_cv_sets.mat');        
    elseif kfold == 5
        load('Prediction_Histology_Result_5-fold_multiclass_healthy_control.mat');        
    end
    
    load('permutation_test_classifier_kfold_multiclass_trial2.mat');
    
end

%% Read demograpic data and set up some parameters and variables
for read_demodata_setup_params = 1
    
    BASEDIR = '/local_raid/seokjun/01_project/04_IntracorticalAnalysis/';
    OUTPATH = [ BASEDIR '03_result/prediction_histology_by_MRI/' ];
    Cases_cont                = [ BASEDIR '/01_analysis/Demographic_data_control.txt' ];
    Cases_cont_dti            = [ BASEDIR '/01_analysis/Demographic_data_control_DTI.txt' ];
    Cases_cont_fMRI           = [ BASEDIR '/01_analysis/Demographic_data_control_fMRI.txt' ];
    Group_cont                = 'control';
    Prefix_cont               = 'TLE';
    
    Cases_pat                 = [ BASEDIR '/01_analysis/Demographic_data_FCD.txt' ];
    Cases_pat_dti             = [ BASEDIR '/01_analysis/Demographic_data_FCD_DTI.txt' ];
    Cases_pat_fMRI            = [ BASEDIR '/01_analysis/Demographic_data_FCD_fMRI.txt' ];
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
    average_surface_dir       = [ BASEDIR '/03_result/average_surfaces/' ];
    
    visualization = 1;
    save_fig = 0;
%     load([ BASEDIR '/01_analysis/pipeline/postprocessing_zscore_visualization/colormap_noel1.mat' ]);

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
for t1_flair_dti = 1
    
    delete(gcp);
    parpool(24);
    
    %% 1. Make the feature vectors, matching the T1, FLAIR and DTI zscore
    for featureset_build = 1
        
        for control_setup = 1
            
            case_num_cont = case_num_cont_org;
            case_num_cont_dti = case_num_cont_dti_org;            
            pl_feature_set_T1_FLAIR_cont = pl_feature_set_T1_FLAIR_cont_org;
            pl_feature_set_DTI_cont = pl_feature_set_DTI_cont_org;            
            mean_z_lesion_T1_FLAIR_cont = mean_z_lesion_T1_FLAIR_cont_org;
            mean_z_lesion_DTI_cont = mean_z_lesion_DTI_cont_org;
            
            %% T1/FLAIR
            excluded_cases = { '306_1', '313_1', '322_1' };
            [C, ia1, ib1] = intersect(case_num_cont, excluded_cases);
            case_num_cont(ia1) = [];
          
            %% DTI
            excluded_cases = { '306_1', '313_1', '322_1' };
            [C, ia1, ib1] = intersect(case_num_cont_dti, excluded_cases);
            case_num_cont_dti(ia1) = [];
           
            %% final combination
            mean_z_perilesion_T1_cont_final    = mean(pl_feature_set_T1_FLAIR_cont(6).data, 4)'; mean_z_perilesion_T1_cont_final = mean(mean_z_perilesion_T1_cont_final(:, 1:3), 2)';
            mean_z_perilesion_FLAIR_cont_final = mean(squeeze(mean(pl_feature_set_T1_FLAIR_cont(9).data(:, :, 2:5, :), 3)), 3)'; mean_z_perilesion_FLAIR_cont_final = mean(mean_z_perilesion_FLAIR_cont_final(:, 1:6), 2)';
            mean_z_perilesion_DTI_cont_final   = mean(squeeze(mean(pl_feature_set_DTI_cont(1).data(:, :, 2:5, :), 3)), 3)'; mean_z_perilesion_DTI_cont_final = mean(mean_z_perilesion_DTI_cont_final(:, 1:5), 2)';
            
            mean_z_lesion_cont      = [ squeeze(mean(mean_z_lesion_T1_FLAIR_cont, 2)); squeeze(mean(mean_z_lesion_DTI_cont, 2)); mean_z_perilesion_T1_cont_final; mean_z_perilesion_FLAIR_cont_final; mean_z_perilesion_DTI_cont_final;  ];
            mean_z_lesion_cont_temp = [ squeeze(mean(mean_z_lesion_T1_FLAIR_cont, 2)); squeeze(mean(mean_z_lesion_DTI_cont, 2)); ];
%             mean_z_lesion_cont_temp(79, 6) = 0.7273;
            
            peri_lesional_dist      =  [ 2:2:4 ];          
            mean_z_lesion_cont_temp2 = [ squeeze(mean(mean_z_lesion_T1_FLAIR_cont, 2)); squeeze(mean(mean_z_lesion_DTI_cont, 2)); 
                                         mean(pl_feature_set_T1_FLAIR_cont(6).data(peri_lesional_dist, :, :, :), 4);                                                                        %% T1 cortical thickness
                                         mean(squeeze(mean(pl_feature_set_T1_FLAIR_cont(9).data(peri_lesional_dist, :, 2:5, :), 3)), 3)                                                     %% FLAIR cortical intentisy
                                         mean(squeeze(mean(pl_feature_set_DTI_cont(1).data(peri_lesional_dist, :, 2:5, :), 3)), 3);                                                         %% DTI cortical FA                                       
                                       ];
            
        end

        for patient_setup = 1
            
            [C, ia1, ib1] = intersect(case_num_pat, case_num_pat_dti);
            
            FCD_IIa_org = strncmp(histo_type, 'FCDIIa', 6)';
            FCD_IIb_org = strncmp(histo_type, 'FCDIIb', 6)';
            FCD_IIa_dti = strncmp(histo_type_dti, 'FCDIIa', 6)';
            FCD_IIb_dti = strncmp(histo_type_dti, 'FCDIIb', 6)';
                       
            mean_z_perilesion_T1_final    = zeros(1, length(case_num_pat));
            mean_z_perilesion_FLAIR_final = zeros(6, length(case_num_pat));
            mean_z_perilesion_DTI_final   = zeros(5, length(case_num_pat_dti));           
            
            mean_z_perilesion_T1_final(FCD_IIa_org) = mean_z_perilesion_T1.mean_z_score_typeIIa';
            mean_z_perilesion_T1_final(FCD_IIb_org) = mean_z_perilesion_T1.mean_z_score_typeIIb';
            mean_z_perilesion_FLAIR_final(:, FCD_IIa_org) = mean_z_perilesion_FLAIR.mean_z_score_typeIIa';
            mean_z_perilesion_FLAIR_final(:, FCD_IIb_org) = mean_z_perilesion_FLAIR.mean_z_score_typeIIb';
            mean_z_perilesion_DTI_final(:, FCD_IIa_dti) = mean_z_perilesion_DTI.mean_z_score_typeIIa';
            mean_z_perilesion_DTI_final(:, FCD_IIb_dti) = mean_z_perilesion_DTI.mean_z_score_typeIIb';          
            
            mean_z_perilesion_T1_final = mean_z_perilesion_T1_final(1, ia1);
            mean_z_perilesion_FLAIR_final = mean_z_perilesion_FLAIR_final(:, ia1);
            
            mean_z_lesion_FCD      = [ mean_z_lesion_T1_FLAIR(:, ia1); mean_z_lesion_DTI; mean(mean_z_perilesion_T1_final, 1); mean(mean_z_perilesion_FLAIR_final, 1); mean(mean_z_perilesion_DTI_final, 1);  ];
            mean_z_lesion_FCD_temp = [ mean_z_lesion_T1_FLAIR(:, ia1); mean_z_lesion_DTI; ];
            
            peri_lesional_dist      =  [ 2:2:4 ];           
            mean_z_lesion_FCD_temp2 = [ mean_z_lesion_T1_FLAIR(:, ia1); mean_z_lesion_DTI;
                                    pl_feature_set_T1_FLAIR(6).data(peri_lesional_dist, ia1);                                                                                          %% T1 cortical thickness
                                    mean(pl_feature_set_T1_FLAIR(9).data(peri_lesional_dist, ia1, 2:5), 3);                                                                            %% FLAIR cortical intentisy
                                    mean(pl_feature_set_DTI(1).data(peri_lesional_dist, :, [3 5]), 3);                                                                                 %% DTI cortical FA                                   
                                  ];
            
        end
        
        mean_z_lesion        = [ mean_z_lesion_cont mean_z_lesion_FCD ];
        mean_z_lesion_temp   = [ mean_z_lesion_cont_temp mean_z_lesion_FCD_temp ];
        mean_z_lesion_temp2  = [ mean_z_lesion_cont_temp2 mean_z_lesion_FCD_temp2 ];
        
        group_class = cell(size(mean_z_lesion_cont, 2)+size(mean_z_lesion_FCD, 2), 1);
        group_class(1:size(mean_z_lesion_cont, 2)) = { 'control' };
        group_class(size(mean_z_lesion_cont, 2)+1:end) = histo_type_dti;
        
    end
    
    %% 2. Predefine k-fold cross-validation sets
    iteration_num = 100; 
    kfold = 10;
    
    for i = 1 : iteration_num
        
        Indices = crossvalind('kfold', group_class, kfold);
        kfold_set(i, :) = Indices';        
        
    end
    
    bagging_num = 500;

    %% 3. Use unimodal global features based on pre-selection
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
        %% 82 83 = Thickness of T1 perilesion    -> 82 83    p24
        %% 84 85 = Intensity of FLAIR perilesion -> 84 85    p25
        %% 86 87 = FA of DTI perilesion          -> 86 87    p26
        
        %% Using T1 global
        for T1_glboal_features = 1
            
            % one vs. one 
            warning off;
            classified_group = group_class;
            SVMstruct_set = cell(length(classified_group), 1);
            discriminative_features_set = cell(length(group_class), 1);      
            discriminative_features = { [2:5], [7:9 16], [12:15], 17, 18, 19, [20:22] };
            opts = statset();      
            
            mean_z_lesion_final = [];
            for i = 1 : length(discriminative_features)
                mean_z_lesion_final = [ mean_z_lesion_final; mean(mean_z_lesion_temp(discriminative_features{i}, :), 1) ];
            end
                        
            keepout_vec = logical(zeros(size(mean_z_lesion_final, 1), 1));
            keepout_vec([  ]) = 1;
            
            CONTROL    = strcmp(group_class, 'control');
            FCDTypeIIA = strcmp(group_class, 'FCDIIa');
            FCDTypeIIB = strcmp(group_class, 'FCDIIb');
            
            prediction_set = zeros(length(group_class), iteration_num);
            prediction_real_set = zeros(length(group_class), iteration_num);
            fs_set = zeros(kfold, size(mean_z_lesion_final, 1), iteration_num);
            
            for i = 1 : iteration_num
                
                Indices = kfold_set(i, :)';
                prediction_curr = zeros(length(group_class), 1);
                prediction_curr_real = zeros(length(group_class), 1);
                
                fs_set_temp = cell(kfold, 1);
                history_set_temp = cell(kfold, 1);
                p_val_temp = cell(kfold, 1);
                
                parfor c = 1 : kfold
                    
                    test_cases = Indices==c;                    
                    
                    CV_group = group_class;
                    CV_group(test_cases) = [];
                    
                    mean_z_lesion_temp_lvo = mean_z_lesion_final;
                    mean_z_lesion_temp_lvo(:, test_cases) = [];
                    
                    p = cvpartition(CV_group,'resubstitution');
                    [fs,history] = sequentialfs(@fs_func_multiclass_svm, mean_z_lesion_temp_lvo', CV_group,'cv', p, 'direction', 'forward', 'keepout', keepout_vec', 'options',opts);
                    fs_set_temp{c} = fs;
                    
                end
                
                fs_set_temp = cell2mat(fs_set_temp);                
                
                for c = 1 : kfold
                    
                    test_cases  = Indices==c;
                    train_cases = Indices~=c;
                    
                    train_control = find(strcmp(group_class, 'control').*train_cases);
                    train_groupA = find(strcmp(group_class, 'FCDIIa').*train_cases);
                    train_groupB = find(strcmp(group_class, 'FCDIIb').*train_cases);
                    
                    num_of_control_in_train = length(train_control);
                    num_of_FCDIIA_in_train = length(train_groupA);
                    num_of_FCDIIB_in_train = length(train_groupB);
                    num_of_each_group = min( [ num_of_control_in_train num_of_FCDIIA_in_train num_of_FCDIIB_in_train] );
                    
                    CV_group = group_class;
                    CV_group(test_cases) = [];
                    
                    mean_z_lesion_temp_lvo = mean_z_lesion_final;
                    mean_z_lesion_temp_lvo(:, test_cases) = [];
                    
                    test_real  = find(test_cases);
                    CV_group_test = group_class(test_real);
                    
                    fs = fs_set_temp(c, :);
                    
                    if(bagging_num==0)
                        
                        classified_group = cell(length(group_class), 1);
                        train_real = find(train_cases);
                        test_real  = find(test_cases);
                        
                        CV_group_train = group_class(train_real);
                        CV_group_test = group_class(test_real);
                        
                        training_data = mean_z_lesion_final(fs, train_real);
                        test_data     = mean_z_lesion_final(fs, test_real);
                        
                        t = templateSVM('KernelFunction', 'polynomial');
                        classOrder = unique(group_class);
                        
                        Mdl = fitcecoc(training_data',CV_group_train, 'Learners', t, 'ClassNames', classOrder);
                        newClass = predict(Mdl, test_data');
                        
                        ia = [];
                        ia(strcmp(newClass, 'control'), 1) = 1;
                        ia(strcmp(newClass, 'FCDIIa'), 1) = 2;
                        ia(strcmp(newClass, 'FCDIIb'), 1) = 3;
                        
                        prediction_curr(test_cases) = strcmp(CV_group_test, newClass);
                        prediction_curr_real(test_cases) = ia;
                        
                    elseif(bagging_num>0)
                        
                        test_bagging = cell(1, bagging_num);
                        CV_group_test = group_class(test_real);
                        
                        parfor j = 1 : bagging_num
                            
                            train_control_real = sort(train_control(randperm(num_of_control_in_train, num_of_each_group)));
                            train_groupA_real  = sort(train_groupA(randperm(num_of_FCDIIA_in_train, num_of_each_group)));
                            train_groupB_real  = sort(train_groupB(randperm(num_of_FCDIIB_in_train, num_of_each_group)));
                            
                            train_real = sort([ train_control_real; train_groupA_real; train_groupB_real]);
                            test_real  = find(test_cases);
                            
                            CV_group_train = group_class(train_real);                            
                            
                            training_data = mean_z_lesion_final(fs, train_real);
                            test_data     = mean_z_lesion_final(fs, test_real);
                            
                            t = templateSVM('KernelFunction', 'polynomial');
                            classOrder = unique(group_class);
                            
                            Mdl = fitcecoc(training_data',CV_group_train, 'Learners', t, 'ClassNames', classOrder);
                            newClass = predict(Mdl, test_data');                            
                            
                            ib = [];
                            ib(strcmp(newClass, 'control'), 1) = 1;
                            ib(strcmp(newClass, 'FCDIIa'), 1) = 2;
                            ib(strcmp(newClass, 'FCDIIb'), 1) = 3;
                            
                            test_bagging{j} = ib;
                            
                        end
                        
                        ia = [];
                        ia(strcmp(CV_group_test, 'control'), 1) = 1;
                        ia(strcmp(CV_group_test, 'FCDIIa'), 1) = 2;
                        ia(strcmp(CV_group_test, 'FCDIIb'), 1) = 3;
                        
                        prediction_curr(test_cases) = mode(double(cell2mat(test_bagging)), 2) == ia;
                        prediction_curr_real(test_cases) = mode(double(cell2mat(test_bagging)), 2);
                        
                    end
                    
                end
                    
                disp([ 'iter ', num2str(i), ' : ', num2str(sum(prediction_curr)/length(group_class))]);
                prediction_set(:, i) = prediction_curr;
                prediction_real_set(:, i) = prediction_curr_real; % 1: control, 2: Type IIA, 3: Type IIB
                fs_set(:, :, i) = fs_set_temp;
                
            end
            
            T1_prediction_set = prediction_set;
            T1_prediction_real_set = prediction_real_set;
            T1_fs_set = fs_set;
            
        end
        
        %% Using FLAIR global
        for FLAIR_glboal_features = 1
                        
            % one vs. one 
            warning off;
            classified_group = group_class;
            SVMstruct_set = cell(length(classified_group), 1);
            discriminative_features_set = cell(length(group_class), 1);      
            discriminative_features = { [30:33], [35:37 44], [40:43], [45:47] };
            opts = statset();      
            
            mean_z_lesion_final = [];
            for i = 1 : length(discriminative_features)
                mean_z_lesion_final = [ mean_z_lesion_final; mean(mean_z_lesion_temp(discriminative_features{i}, :), 1) ];
            end
                        
            keepout_vec = logical(zeros(size(mean_z_lesion_final, 1), 1));
            keepout_vec([  ]) = 1;
            
            CONTROL    = strcmp(group_class, 'control');
            FCDTypeIIA = strcmp(group_class, 'FCDIIa');
            FCDTypeIIB = strcmp(group_class, 'FCDIIb');
            
            prediction_set = zeros(length(group_class), iteration_num);
            prediction_real_set = zeros(length(group_class), iteration_num);
            fs_set = zeros(kfold, size(mean_z_lesion_final, 1), iteration_num);            
            
            for i = 1 : iteration_num
                
                Indices = kfold_set(i, :)';
                prediction_curr = zeros(length(group_class), 1);
                prediction_curr_real = zeros(length(group_class), 1);
                
                fs_set_temp = cell(kfold, 1);
                history_set_temp = cell(kfold, 1);
                p_val_temp = cell(kfold, 1);
                
                parfor c = 1 : kfold
                    
                    test_cases = Indices==c;                    
                    
                    CV_group = group_class;
                    CV_group(test_cases) = [];
                    
                    mean_z_lesion_temp_lvo = mean_z_lesion_final;
                    mean_z_lesion_temp_lvo(:, test_cases) = [];
                    
                    p = cvpartition(CV_group,'resubstitution');
                    [fs,history] = sequentialfs(@fs_func_multiclass_svm, mean_z_lesion_temp_lvo', CV_group,'cv', p, 'direction', 'forward', 'keepout', keepout_vec', 'options',opts);
                    fs_set_temp{c} = fs;
                    
                end
                
                fs_set_temp = cell2mat(fs_set_temp);                
                
                for c = 1 : kfold
                    
                    test_cases  = Indices==c;
                    train_cases = Indices~=c;
                    
                    train_control = find(strcmp(group_class, 'control').*train_cases);
                    train_groupA = find(strcmp(group_class, 'FCDIIa').*train_cases);
                    train_groupB = find(strcmp(group_class, 'FCDIIb').*train_cases);
                    
                    num_of_control_in_train = length(train_control);
                    num_of_FCDIIA_in_train = length(train_groupA);
                    num_of_FCDIIB_in_train = length(train_groupB);
                    num_of_each_group = min( [ num_of_control_in_train num_of_FCDIIA_in_train num_of_FCDIIB_in_train] );
                    
                    CV_group = group_class;
                    CV_group(test_cases) = [];
                    
                    mean_z_lesion_temp_lvo = mean_z_lesion_final;
                    mean_z_lesion_temp_lvo(:, test_cases) = [];
                    
                    test_real  = find(test_cases);
                    CV_group_test = group_class(test_real);
                    
                    fs = fs_set_temp(c, :);
                    
                    if(bagging_num==0)
                        
                        classified_group = cell(length(group_class), 1);
                        train_real = find(train_cases);
                        test_real  = find(test_cases);
                        
                        CV_group_train = group_class(train_real);
                        CV_group_test = group_class(test_real);
                        
                        training_data = mean_z_lesion_final(fs, train_real);
                        test_data     = mean_z_lesion_final(fs, test_real);
                        
                        t = templateSVM('KernelFunction', 'polynomial', 'Standardize', 1, 'KernelScale', 'auto');
                        classOrder = unique(group_class);
                        
                        Mdl = fitcecoc(training_data',CV_group_train, 'Learners', t, 'ClassNames', classOrder);
                        newClass = predict(Mdl, test_data');
                        
                        ia = [];
                        ia(strcmp(newClass, 'control'), 1) = 1;
                        ia(strcmp(newClass, 'FCDIIa'), 1) = 2;
                        ia(strcmp(newClass, 'FCDIIb'), 1) = 3;
                        
                        prediction_curr(test_cases) = strcmp(CV_group_test, newClass);
                        prediction_curr_real(test_cases) = ia;
                        
                    elseif(bagging_num>0)
                        
                        test_bagging = cell(1, bagging_num);
                        CV_group_test = group_class(test_real);
                        
                        parfor j = 1 : bagging_num
                            
                            train_control_real = sort(train_control(randperm(num_of_control_in_train, num_of_each_group)));
                            train_groupA_real  = sort(train_groupA(randperm(num_of_FCDIIA_in_train, num_of_each_group)));
                            train_groupB_real  = sort(train_groupB(randperm(num_of_FCDIIB_in_train, num_of_each_group)));
                            
                            train_real = sort([ train_control_real; train_groupA_real; train_groupB_real]);
                            test_real  = find(test_cases);
                            
                            CV_group_train = group_class(train_real);                            
                            
                            training_data = mean_z_lesion_final(fs, train_real);
                            test_data     = mean_z_lesion_final(fs, test_real);
                            
                            t = templateSVM('KernelFunction', 'polynomial', 'Standardize', 1, 'KernelScale', 'auto');
                            classOrder = unique(group_class);
                            
                            Mdl = fitcecoc(training_data',CV_group_train, 'Learners', t, 'ClassNames', classOrder);
                            newClass = predict(Mdl, test_data');                            
                            
                            ib = [];
                            ib(strcmp(newClass, 'control'), 1) = 1;
                            ib(strcmp(newClass, 'FCDIIa'), 1) = 2;
                            ib(strcmp(newClass, 'FCDIIb'), 1) = 3;
                            
                            test_bagging{j} = ib;
                            
                        end
                        
                        ia = [];
                        ia(strcmp(CV_group_test, 'control'), 1) = 1;
                        ia(strcmp(CV_group_test, 'FCDIIa'), 1) = 2;
                        ia(strcmp(CV_group_test, 'FCDIIb'), 1) = 3;
                        
                        prediction_curr(test_cases) = mode(double(cell2mat(test_bagging)), 2) == ia;
                        prediction_curr_real(test_cases) = mode(double(cell2mat(test_bagging)), 2);
                        
                    end
                    
                end
                    
                disp([ 'iter ', num2str(i), ' : ', num2str(sum(prediction_curr)/length(group_class))]);
                prediction_set(:, i) = prediction_curr;
                prediction_real_set(:, i) = prediction_curr_real; % 1: control, 2: Type IIA, 3: Type IIB
                fs_set(:, :, i) = fs_set_temp;
                
            end
            
            FLAIR_prediction_set = prediction_set;
            FLAIR_prediction_real_set = prediction_real_set;
            FLAIR_fs_set = fs_set;
            
        end
        
        %% Using DTI global
        for DTI_glboal_features = 1
                      
            % one vs. one 
            warning off;
            classified_group = group_class;
            SVMstruct_set = cell(length(classified_group), 1);
            discriminative_features_set = cell(length(group_class), 1);      
            discriminative_features = { [56 58], [60 62], [68 70], [72 74] };
            opts = statset();      
            
            mean_z_lesion_final = [];
            for i = 1 : length(discriminative_features)
                mean_z_lesion_final = [ mean_z_lesion_final; mean(mean_z_lesion_temp(discriminative_features{i}, :), 1) ];
            end
                        
            keepout_vec = logical(zeros(size(mean_z_lesion_final, 1), 1));
            keepout_vec([  ]) = 1;
            
            CONTROL    = strcmp(group_class, 'control');
            FCDTypeIIA = strcmp(group_class, 'FCDIIa');
            FCDTypeIIB = strcmp(group_class, 'FCDIIb');
            
            prediction_set = zeros(length(group_class), iteration_num);
            prediction_real_set = zeros(length(group_class), iteration_num);
            fs_set = zeros(kfold, size(mean_z_lesion_final, 1), iteration_num);
            
            for i = 1 : iteration_num
                
                Indices = kfold_set(i, :)';
                prediction_curr = zeros(length(group_class), 1);
                prediction_curr_real = zeros(length(group_class), 1);
                
                fs_set_temp = cell(kfold, 1);
                history_set_temp = cell(kfold, 1);
                p_val_temp = cell(kfold, 1);
                
                parfor c = 1 : kfold
                    
                    test_cases = Indices==c;                    
                    
                    CV_group = group_class;
                    CV_group(test_cases) = [];
                    
                    mean_z_lesion_temp_lvo = mean_z_lesion_final;
                    mean_z_lesion_temp_lvo(:, test_cases) = [];
                    
                    p = cvpartition(CV_group,'resubstitution');
                    [fs,history] = sequentialfs(@fs_func_multiclass_svm, mean_z_lesion_temp_lvo', CV_group,'cv', p, 'direction', 'forward', 'keepout', keepout_vec', 'options',opts);
                    fs_set_temp{c} = fs;
                    
                end
                
                fs_set_temp = cell2mat(fs_set_temp);                
                
                for c = 1 : kfold
                    
                    test_cases  = Indices==c;
                    train_cases = Indices~=c;
                    
                    train_control = find(strcmp(group_class, 'control').*train_cases);
                    train_groupA = find(strcmp(group_class, 'FCDIIa').*train_cases);
                    train_groupB = find(strcmp(group_class, 'FCDIIb').*train_cases);
                    
                    num_of_control_in_train = length(train_control);
                    num_of_FCDIIA_in_train = length(train_groupA);
                    num_of_FCDIIB_in_train = length(train_groupB);
                    num_of_each_group = min( [ num_of_control_in_train num_of_FCDIIA_in_train num_of_FCDIIB_in_train] );
                    
                    CV_group = group_class;
                    CV_group(test_cases) = [];
                    
                    mean_z_lesion_temp_lvo = mean_z_lesion_final;
                    mean_z_lesion_temp_lvo(:, test_cases) = [];
                    
                    test_real  = find(test_cases);
                    CV_group_test = group_class(test_real);
                    
                    fs = fs_set_temp(c, :);
                    
                    if(bagging_num==0)
                        
                        classified_group = cell(length(group_class), 1);
                        train_real = find(train_cases);
                        test_real  = find(test_cases);
                        
                        CV_group_train = group_class(train_real);
                        CV_group_test = group_class(test_real);
                        
                        training_data = mean_z_lesion_final(fs, train_real);
                        test_data     = mean_z_lesion_final(fs, test_real);
                        
                        t = templateSVM('KernelFunction', 'polynomial');
                        classOrder = unique(group_class);
                        
                        Mdl = fitcecoc(training_data',CV_group_train, 'Learners', t, 'ClassNames', classOrder);
                        newClass = predict(Mdl, test_data');
                        
                        ia = [];
                        ia(strcmp(newClass, 'control'), 1) = 1;
                        ia(strcmp(newClass, 'FCDIIa'), 1) = 2;
                        ia(strcmp(newClass, 'FCDIIb'), 1) = 3;
                        
                        prediction_curr(test_cases) = strcmp(CV_group_test, newClass);
                        prediction_curr_real(test_cases) = ia;
                        
                    elseif(bagging_num>0)
                        
                        test_bagging = cell(1, bagging_num);
                        CV_group_test = group_class(test_real);
                        
                        parfor j = 1 : bagging_num
                            
                            train_control_real = sort(train_control(randperm(num_of_control_in_train, num_of_each_group)));
                            train_groupA_real  = sort(train_groupA(randperm(num_of_FCDIIA_in_train, num_of_each_group)));
                            train_groupB_real  = sort(train_groupB(randperm(num_of_FCDIIB_in_train, num_of_each_group)));
                            
                            train_real = sort([ train_control_real; train_groupA_real; train_groupB_real]);
                            test_real  = find(test_cases);
                            
                            CV_group_train = group_class(train_real);                            
                            
                            training_data = mean_z_lesion_final(fs, train_real);
                            test_data     = mean_z_lesion_final(fs, test_real);
                            
                            t = templateSVM('KernelFunction', 'polynomial');
                            classOrder = unique(group_class);
                            
                            Mdl = fitcecoc(training_data',CV_group_train, 'Learners', t, 'ClassNames', classOrder);
                            newClass = predict(Mdl, test_data');                            
                            
                            ib = [];
                            ib(strcmp(newClass, 'control'), 1) = 1;
                            ib(strcmp(newClass, 'FCDIIa'), 1) = 2;
                            ib(strcmp(newClass, 'FCDIIb'), 1) = 3;
                            
                            test_bagging{j} = ib;
                            
                        end
                        
                        ia = [];
                        ia(strcmp(CV_group_test, 'control'), 1) = 1;
                        ia(strcmp(CV_group_test, 'FCDIIa'), 1) = 2;
                        ia(strcmp(CV_group_test, 'FCDIIb'), 1) = 3;
                        
                        prediction_curr(test_cases) = mode(double(cell2mat(test_bagging)), 2) == ia;
                        prediction_curr_real(test_cases) = mode(double(cell2mat(test_bagging)), 2);
                        
                    end
                    
                end
                    
                disp([ 'iter ', num2str(i), ' : ', num2str(sum(prediction_curr)/length(group_class))]);
                prediction_set(:, i) = prediction_curr;
                prediction_real_set(:, i) = prediction_curr_real; % 1: control, 2: Type IIA, 3: Type IIB
                fs_set(:, :, i) = fs_set_temp;
                
            end
            
            DTI_prediction_set = prediction_set;
            DTI_prediction_real_set = prediction_real_set;
            DTI_fs_set = fs_set;
            
        end
        
    end
    
    %% 4. Use multimodal global features based on pre-selection
    for prediction_using_multimodal_global_features = 1
        
        %% Using multimodal global features
        for multimodal_glboal_features = 1
                        
            % one vs. one 
            warning off;
            classified_group = group_class;
            SVMstruct_set = cell(length(classified_group), 1);
            discriminative_features_set = cell(length(group_class), 1);      
            discriminative_features = { [2:5], [7:9 16], [12:15], 17, 18, 19, [20:22], [30:33], [35:37 44], [40:43], [45:47], [56 58], [60 62], [68 70], [72 74] };
            opts = statset();      
            
            mean_z_lesion_final = [];
            for i = 1 : length(discriminative_features)
                mean_z_lesion_final = [ mean_z_lesion_final; mean(mean_z_lesion_temp(discriminative_features{i}, :), 1) ];
            end
                        
            keepout_vec = logical(zeros(size(mean_z_lesion_final, 1), 1));
            keepout_vec([  ]) = 1;
            
            CONTROL    = strcmp(group_class, 'control');
            FCDTypeIIA = strcmp(group_class, 'FCDIIa');
            FCDTypeIIB = strcmp(group_class, 'FCDIIb');
            
            prediction_set = zeros(length(group_class), iteration_num);
            prediction_real_set = zeros(length(group_class), iteration_num);
            fs_set = zeros(kfold, size(mean_z_lesion_final, 1), iteration_num);
            
            for i = 1 : iteration_num
                
                Indices = kfold_set(i, :)';
                prediction_curr = zeros(length(group_class), 1);
                prediction_curr_real = zeros(length(group_class), 1);
                
                fs_set_temp = cell(kfold, 1);
                history_set_temp = cell(kfold, 1);
                p_val_temp = cell(kfold, 1);
                
                parfor c = 1 : kfold
                    
                    test_cases = Indices==c;                    
                    
                    CV_group = group_class;
                    CV_group(test_cases) = [];
                    
                    mean_z_lesion_temp_lvo = mean_z_lesion_final;
                    mean_z_lesion_temp_lvo(:, test_cases) = [];
                    
                    p = cvpartition(CV_group,'resubstitution');
                    [fs,history] = sequentialfs(@fs_func_multiclass_svm, mean_z_lesion_temp_lvo', CV_group,'cv', p, 'direction', 'forward', 'keepout', keepout_vec', 'options',opts);
                    fs_set_temp{c} = fs;
                    
                end
                
                fs_set_temp = cell2mat(fs_set_temp);                
                
                for c = 1 : kfold
                    
                    test_cases  = Indices==c;
                    train_cases = Indices~=c;
                    
                    train_control = find(strcmp(group_class, 'control').*train_cases);
                    train_groupA = find(strcmp(group_class, 'FCDIIa').*train_cases);
                    train_groupB = find(strcmp(group_class, 'FCDIIb').*train_cases);
                    
                    num_of_control_in_train = length(train_control);
                    num_of_FCDIIA_in_train = length(train_groupA);
                    num_of_FCDIIB_in_train = length(train_groupB);
                    num_of_each_group = min( [ num_of_control_in_train num_of_FCDIIA_in_train num_of_FCDIIB_in_train] );
                    
                    CV_group = group_class;
                    CV_group(test_cases) = [];
                    
                    mean_z_lesion_temp_lvo = mean_z_lesion_final;
                    mean_z_lesion_temp_lvo(:, test_cases) = [];
                    
                    test_real  = find(test_cases);
                    CV_group_test = group_class(test_real);
                    
                    fs = fs_set_temp(c, :);
                    
                    if(bagging_num==0)
                        
                        classified_group = cell(length(group_class), 1);
                        train_real = find(train_cases);
                        test_real  = find(test_cases);
                        
                        CV_group_train = group_class(train_real);
                        CV_group_test = group_class(test_real);
                        
                        training_data = mean_z_lesion_final(fs, train_real);
                        test_data     = mean_z_lesion_final(fs, test_real);
                        
                        t = templateSVM('KernelFunction', 'polynomial');
                        classOrder = unique(group_class);
                        
                        Mdl = fitcecoc(training_data',CV_group_train, 'Learners', t, 'ClassNames', classOrder);
                        newClass = predict(Mdl, test_data');
                        
                        ia = [];
                        ia(strcmp(newClass, 'control'), 1) = 1;
                        ia(strcmp(newClass, 'FCDIIa'), 1) = 2;
                        ia(strcmp(newClass, 'FCDIIb'), 1) = 3;
                        
                        prediction_curr(test_cases) = strcmp(CV_group_test, newClass);
                        prediction_curr_real(test_cases) = ia;
                        
                    elseif(bagging_num>0)
                        
                        test_bagging = cell(1, bagging_num);
                        CV_group_test = group_class(test_real);
                        
                        parfor j = 1 : bagging_num
                            
                            train_control_real = sort(train_control(randperm(num_of_control_in_train, num_of_each_group)));
                            train_groupA_real  = sort(train_groupA(randperm(num_of_FCDIIA_in_train, num_of_each_group)));
                            train_groupB_real  = sort(train_groupB(randperm(num_of_FCDIIB_in_train, num_of_each_group)));
                            
                            train_real = sort([ train_control_real; train_groupA_real; train_groupB_real]);
                            test_real  = find(test_cases);
                            
                            CV_group_train = group_class(train_real);                            
                            
                            training_data = mean_z_lesion_final(fs, train_real);
                            test_data     = mean_z_lesion_final(fs, test_real);
                            
                            t = templateSVM('KernelFunction', 'polynomial');
                            classOrder = unique(group_class);
                            
                            Mdl = fitcecoc(training_data',CV_group_train, 'Learners', t, 'ClassNames', classOrder);
                            newClass = predict(Mdl, test_data');                            
                            
                            ib = [];
                            ib(strcmp(newClass, 'control'), 1) = 1;
                            ib(strcmp(newClass, 'FCDIIa'), 1) = 2;
                            ib(strcmp(newClass, 'FCDIIb'), 1) = 3;
                            
                            test_bagging{j} = ib;
                            
                        end
                        
                        ia = [];
                        ia(strcmp(CV_group_test, 'control'), 1) = 1;
                        ia(strcmp(CV_group_test, 'FCDIIa'), 1) = 2;
                        ia(strcmp(CV_group_test, 'FCDIIb'), 1) = 3;
                        
                        prediction_curr(test_cases) = mode(double(cell2mat(test_bagging)), 2) == ia;
                        prediction_curr_real(test_cases) = mode(double(cell2mat(test_bagging)), 2);
                        
                    end
                    
                end
                    
                disp([ 'iter ', num2str(i), ' : ', num2str(sum(prediction_curr)/length(group_class))]);
                prediction_set(:, i) = prediction_curr;
                prediction_real_set(:, i) = prediction_curr_real; % 1: control, 2: Type IIA, 3: Type IIB
                fs_set(:, :, i) = fs_set_temp;
                
            end
            
            mMGlobal_prediction_set = prediction_set;
            mMGlobal_prediction_real_set = prediction_real_set;
            mMGlobal_fs_set = fs_set;
            
        end 
        
    end
    
    %% 5. Use depth multi-modality features based on pre-selection
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
        %% 78 79 = Thickness of T1 perilesion    -> 78 79    p24
        %% 80 81 = Intensity of FLAIR perilesion -> 80 81    p25
        %% 82 83 = FA of DTI perilesion          -> 82 83    p26
        
        %% Using multimodal + multisurface 
        for all_features_multimodal_multisurface = 1

            % one vs. one 
            warning off;
            classified_group = group_class;
            SVMstruct_set = cell(length(classified_group), 1);
            discriminative_features_set = cell(length(group_class), 1);            
            opts = statset();            
            mean_z_lesion_final = mean_z_lesion_temp;
            keepout_vec = logical(zeros(size(mean_z_lesion_final, 1), 1));
            keepout_vec([ 1 6 10 11 23:28 29 34 38 39 48:53 54 55 57 59 61 63:65 66 67 69 71 73 75:77]) = 1;            
            
            CONTROL    = strcmp(group_class, 'control');
            FCDTypeIIA = strcmp(group_class, 'FCDIIa');
            FCDTypeIIB = strcmp(group_class, 'FCDIIb');
            
            prediction_set = zeros(length(group_class), iteration_num);
            prediction_real_set = zeros(length(group_class), iteration_num);
            fs_set = zeros(kfold, size(mean_z_lesion_final, 1), iteration_num);
            
            for i = 1 : iteration_num
                
                Indices = kfold_set(i, :)';
                prediction_curr = zeros(length(group_class), 1);
                prediction_curr_real = zeros(length(group_class), 1);
                
                fs_set_temp = cell(kfold, 1);
                history_set_temp = cell(kfold, 1);
                p_val_temp = cell(kfold, 1);
                
                parfor c = 1 : kfold
                    
                    test_cases = Indices==c;
                    
                    CV_group = group_class;
                    CV_group(test_cases) = [];
                    
                    mean_z_lesion_temp_lvo = mean_z_lesion_final;
                    mean_z_lesion_temp_lvo(:, test_cases) = [];
                    
                    p = cvpartition(CV_group,'resubstitution');
                    [fs,history] = sequentialfs(@fs_func_multiclass_svm, mean_z_lesion_temp_lvo', CV_group,'cv', p, 'direction', 'forward', 'keepout', keepout_vec', 'options',opts);
                    fs_set_temp{c} = fs;
                    
                end
                
                fs_set_temp = cell2mat(fs_set_temp);
                
                for c = 1 : kfold
                    
                    test_cases  = Indices==c;
                    train_cases = Indices~=c;
                    
                    train_control = find(strcmp(group_class, 'control').*train_cases);
                    train_groupA = find(strcmp(group_class, 'FCDIIa').*train_cases);
                    train_groupB = find(strcmp(group_class, 'FCDIIb').*train_cases);
                    
                    num_of_control_in_train = length(train_control);
                    num_of_FCDIIA_in_train = length(train_groupA);
                    num_of_FCDIIB_in_train = length(train_groupB);
                    num_of_each_group = min( [ num_of_control_in_train num_of_FCDIIA_in_train num_of_FCDIIB_in_train] );
                    
                    CV_group = group_class;
                    CV_group(test_cases) = [];
                    
                    mean_z_lesion_temp_lvo = mean_z_lesion_final;
                    mean_z_lesion_temp_lvo(:, test_cases) = [];
                    
                    test_real  = find(test_cases);
                    CV_group_test = group_class(test_real);
                    
                    fs = fs_set_temp(c, :);
                    
                    if(bagging_num==0)
                        
                        classified_group = cell(length(group_class), 1);
                        train_real = find(train_cases);
                        test_real  = find(test_cases);
                        
                        CV_group_train = group_class(train_real);
                        CV_group_test = group_class(test_real);
                        
                        training_data = mean_z_lesion_final(fs, train_real);
                        test_data     = mean_z_lesion_final(fs, test_real);
                        
                        t = templateSVM('KernelFunction', 'polynomial');
                        classOrder = unique(group_class);
                        
                        Mdl = fitcecoc(training_data',CV_group_train, 'Learners', t, 'ClassNames', classOrder);
                        newClass = predict(Mdl, test_data');
                        
                        ia = [];
                        ia(strcmp(newClass, 'control'), 1) = 1;
                        ia(strcmp(newClass, 'FCDIIa'), 1) = 2;
                        ia(strcmp(newClass, 'FCDIIb'), 1) = 3;
                        
                        prediction_curr(test_cases) = strcmp(CV_group_test, newClass);
                        prediction_curr_real(test_cases) = ia;
                        
                    elseif(bagging_num>0)
                        
                        test_bagging = cell(1, bagging_num);
                        CV_group_test = group_class(test_real);
                        
                        parfor j = 1 : bagging_num
                            
                            train_control_real = sort(train_control(randperm(num_of_control_in_train, num_of_each_group)));
                            train_groupA_real  = sort(train_groupA(randperm(num_of_FCDIIA_in_train, num_of_each_group)));
                            train_groupB_real  = sort(train_groupB(randperm(num_of_FCDIIB_in_train, num_of_each_group)));
                            
                            train_real = sort([ train_control_real; train_groupA_real; train_groupB_real]);
                            test_real  = find(test_cases);
                            
                            CV_group_train = group_class(train_real);
                            
                            training_data = mean_z_lesion_final(fs, train_real);
                            test_data     = mean_z_lesion_final(fs, test_real);
                            
                            t = templateSVM('KernelFunction', 'polynomial');
                            classOrder = unique(group_class);
                            
                            Mdl = fitcecoc(training_data',CV_group_train, 'Learners', t, 'ClassNames', classOrder);
                            newClass = predict(Mdl, test_data');
                            
                            ib = [];
                            ib(strcmp(newClass, 'control'), 1) = 1;
                            ib(strcmp(newClass, 'FCDIIa'), 1) = 2;
                            ib(strcmp(newClass, 'FCDIIb'), 1) = 3;
                            
                            test_bagging{j} = ib;
                            
                        end
                        
                        ia = [];
                        ia(strcmp(CV_group_test, 'control'), 1) = 1;
                        ia(strcmp(CV_group_test, 'FCDIIa'), 1) = 2;
                        ia(strcmp(CV_group_test, 'FCDIIb'), 1) = 3;
                        
                        prediction_curr(test_cases) = mode(double(cell2mat(test_bagging)), 2) == ia;
                        prediction_curr_real(test_cases) = mode(double(cell2mat(test_bagging)), 2);
                        
                    end
                    
                end
                
                disp([ 'iter ', num2str(i), ' : ', num2str(sum(prediction_curr)/length(group_class))]);
                prediction_set(:, i) = prediction_curr;
                prediction_real_set(:, i) = prediction_curr_real; % 1: control, 2: Type IIA, 3: Type IIB
                fs_set(:, :, i) = fs_set_temp;
                
            end
            
            delete(gcp);
            
            mMmS_prediction_set = prediction_set;
            mMmS_prediction_real_set = prediction_real_set;
            mMmS_fs_set = fs_set;
            
        end
        
        %% Using multimodal + multisurface + peri-lesional
        for all_features_multimodal_multisurface_perilesional = 1
            
            % one vs. one 
            warning off;
            classified_group = group_class;
            SVMstruct_set = cell(length(classified_group), 1);
            discriminative_features_set = cell(length(group_class), 1);            
            opts = statset();            
            mean_z_lesion_final = mean_z_lesion_temp2;
            keepout_vec = logical(zeros(size(mean_z_lesion_final, 1), 1));
            keepout_vec([ 1 6 11 23:28 29 34 39 48:53 54 55 57 61 63:65 66 67 69 71 73 75:77]) = 1;            
            
            CONTROL    = strcmp(group_class, 'control');
            FCDTypeIIA = strcmp(group_class, 'FCDIIa');
            FCDTypeIIB = strcmp(group_class, 'FCDIIb');
            
            prediction_set = zeros(length(group_class), iteration_num);
            prediction_real_set = zeros(length(group_class), iteration_num);
            fs_set = zeros(kfold, size(mean_z_lesion_final, 1), iteration_num);
            
            for i = 1 : iteration_num
                
                Indices = kfold_set(i, :)';
                prediction_curr = zeros(length(group_class), 1);
                prediction_curr_real = zeros(length(group_class), 1);
                
                fs_set_temp = cell(kfold, 1);
                history_set_temp = cell(kfold, 1);
                p_val_temp = cell(kfold, 1);
                
                parfor c = 1 : kfold
                    
                    test_cases = Indices==c;
                    
                    CV_group = group_class;
                    CV_group(test_cases) = [];
                    
                    mean_z_lesion_temp_lvo = mean_z_lesion_final;
                    mean_z_lesion_temp_lvo(:, test_cases) = [];
                    
                    p = cvpartition(CV_group,'resubstitution');
                    keepin_vec  = logical(zeros(size(mean_z_lesion_final, 1), 1));
                    keepin_vec([ 2 9 10 16 17 19 20 30 31 32 33 37 38 45 58 59 60  ]) = 1;
                    [fs,history] = sequentialfs(@fs_func_multiclass_svm, mean_z_lesion_temp_lvo', CV_group,'cv', p, 'direction', 'forward', 'keepout', keepout_vec', 'keepin', keepin_vec', 'options',opts);
                    fs_set_temp{c} = fs;
                    
                end
                
                fs_set_temp = cell2mat(fs_set_temp);
                
                for c = 1 : kfold
                    
                    c
                    test_cases  = Indices==c;
                    train_cases = Indices~=c;
                    
                    train_control = find(strcmp(group_class, 'control').*train_cases);
                    train_groupA = find(strcmp(group_class, 'FCDIIa').*train_cases);
                    train_groupB = find(strcmp(group_class, 'FCDIIb').*train_cases);
                    
                    num_of_control_in_train = length(train_control);
                    num_of_FCDIIA_in_train = length(train_groupA);
                    num_of_FCDIIB_in_train = length(train_groupB);
                    num_of_each_group = min( [ num_of_control_in_train num_of_FCDIIA_in_train num_of_FCDIIB_in_train] );
                    
                    CV_group = group_class;
                    CV_group(test_cases) = [];
                    
                    mean_z_lesion_temp_lvo = mean_z_lesion_final;
                    mean_z_lesion_temp_lvo(:, test_cases) = [];
                    
                    test_real  = find(test_cases);
                    CV_group_test = group_class(test_real);
                    
                    fs = fs_set_temp(c, :);                     
                    
                    if(bagging_num==0)
                        
                        classified_group = cell(length(group_class), 1);
                        train_real = find(train_cases);
                        test_real  = find(test_cases);
                        
                        CV_group_train = group_class(train_real);
                        CV_group_test = group_class(test_real);
                        
                        training_data = mean_z_lesion_final(fs, train_real);
                        test_data     = mean_z_lesion_final(fs, test_real);
                        
                        t = templateSVM('KernelFunction', 'polynomial');
                        classOrder = unique(group_class);
                        
                        Mdl = fitcecoc(training_data',CV_group_train, 'Learners', t, 'ClassNames', classOrder);
                        newClass = predict(Mdl, test_data');
                        
                        ia = [];
                        ia(strcmp(newClass, 'control'), 1) = 1;
                        ia(strcmp(newClass, 'FCDIIa'), 1) = 2;
                        ia(strcmp(newClass, 'FCDIIb'), 1) = 3;
                        
                        prediction_curr(test_cases) = strcmp(CV_group_test, newClass);
                        prediction_curr_real(test_cases) = ia;
                        
                    elseif(bagging_num>0)
                        
                        test_bagging = cell(1, bagging_num);
                        CV_group_test = group_class(test_real);
                        
                        parfor j = 1 : bagging_num
                            
                            train_control_real = sort(train_control(randperm(num_of_control_in_train, num_of_each_group)));
                            train_groupA_real  = sort(train_groupA(randperm(num_of_FCDIIA_in_train, num_of_each_group)));
                            train_groupB_real  = sort(train_groupB(randperm(num_of_FCDIIB_in_train, num_of_each_group)));
                            
                            train_real = sort([ train_control_real; train_groupA_real; train_groupB_real]);
                            test_real  = find(test_cases);
                            
                            CV_group_train = group_class(train_real);
                            
                            training_data = mean_z_lesion_final(fs, train_real);
                            test_data     = mean_z_lesion_final(fs, test_real);
                            
                            t = templateSVM('KernelFunction', 'polynomial');
                            classOrder = unique(group_class);
                            
                            Mdl = fitcecoc(training_data',CV_group_train, 'Learners', t, 'ClassNames', classOrder);
                            newClass = predict(Mdl, test_data');
                            
                            ib = [];
                            ib(strcmp(newClass, 'control'), 1) = 1;
                            ib(strcmp(newClass, 'FCDIIa'), 1) = 2;
                            ib(strcmp(newClass, 'FCDIIb'), 1) = 3;
                            
                            test_bagging{j} = ib;
                            
                        end
                        
                        ia = [];
                        ia(strcmp(CV_group_test, 'control'), 1) = 1;
                        ia(strcmp(CV_group_test, 'FCDIIa'), 1) = 2;
                        ia(strcmp(CV_group_test, 'FCDIIb'), 1) = 3;
                        
                        prediction_curr(test_cases) = mode(double(cell2mat(test_bagging)), 2) == ia;
                        prediction_curr_real(test_cases) = mode(double(cell2mat(test_bagging)), 2);
                        
                    end
                    
                end
                
                disp([ 'iter ', num2str(i), ' : ', num2str(sum(prediction_curr)/length(group_class))]);
                prediction_set(:, i) = prediction_curr;
                prediction_real_set(:, i) = prediction_curr_real; % 1: control, 2: Type IIA, 3: Type IIB
                fs_set(:, :, i) = fs_set_temp;
                
            end
            
            %% prediction: mean±SD=0.8532±0.02
            
            delete(gcp);
            
            mMmSPeri_prediction_set = prediction_set;
            mMmSPeri_prediction_real_set = prediction_real_set;
            mMmSPeri_fs_set = fs_set;
            
        end
        
    end
    
    %% 6. Mcnemar test to compare classifier performance (take median performance from each classifier)
    for mcnemar_test = 1
        
        control_idx = find(strcmp(group_class, 'control'));
        FCD_idx = [ (max(control_idx) + 1) : 52 ]';
        T1_median_idx = find(median(mean(T1_prediction_set, 1)) == mean(T1_prediction_set, 1));
        FLAIR_median_idx = find(median(mean(FLAIR_prediction_set, 1)) == mean(FLAIR_prediction_set, 1));
        DTI_median_idx = find(median(mean(DTI_prediction_set, 1)) == mean(DTI_prediction_set, 1));
        fMRI_median_idx = find(median(mean(fMRI_prediction_set, 1)) == mean(fMRI_prediction_set, 1));
        mMGlobal_median_idx = find(abs(median(mean(mMGlobal_prediction_set, 1)) - mean(mMGlobal_prediction_set, 1))<0.01); % mMGlobalf_median_idx = find(abs(mean(mMGlobalf_prediction_set, 1)-0.75)<0.01);
        mMGlobalf_median_idx = find(abs(median(mean(mMGlobalf_prediction_set, 1)) - mean(mMGlobalf_prediction_set, 1))<0.01); % mMGlobalf_median_idx = find(abs(mean(mMGlobalf_prediction_set, 1)-0.75)<0.01);
        mMmS_median_idx = find(median(mean(mMmS_prediction_set, 1)) == mean(mMmS_prediction_set, 1));
        mMmSPeri_median_idx = find(median(mean(mMmSPeri_prediction_set, 1)) == mean(mMmSPeri_prediction_set, 1)); % mMmSPeri_median_idx = find(abs(mean(mMmSPeri_prediction_set, 1)-0.9231)<0.01);
        
        feature1 = 'mMmSPeri';
        feature2 = 'FLAIR';
        eval([ 'feature1_idx = ' feature1 '_median_idx;' ]);
        eval([ 'feature2_idx = ' feature2 '_median_idx;' ]);
        eval([ 'feature1_prediction = ' feature1 '_prediction_set;' ]);
        eval([ 'feature2_prediction = ' feature2 '_prediction_set;' ]);
        
        chi   = zeros(length(feature1_idx), length(feature2_idx));
        p     = zeros(length(feature1_idx), length(feature2_idx));
        alpha = zeros(length(feature1_idx), length(feature2_idx));
        Zb    = zeros(length(feature1_idx), length(feature2_idx));
        pwr   = zeros(length(feature1_idx), length(feature2_idx));
        for i = 1 : length(feature1_idx)
            
            classifier1 = feature1_prediction(:, feature1_idx(i));
            for j = 1 : length(feature2_idx)
                
                classifier2 = feature2_prediction(:, feature2_idx(j));
                
                mcnemar_mat = [ sum(classifier1 == 1 & classifier2 == 1) sum(classifier1 == 1 & classifier2 == 0); ...
                    sum(classifier1 == 0 & classifier2 == 1) sum(classifier1 == 0 & classifier2 == 0) ];
                
                [chi_, p_, alpha_, Zb_, pwr_] = mcnemar(mcnemar_mat);                
                chi(i, j) = chi_; p(i, j) = p_; alpha(i, j) = alpha_; Zb(i, j) = Zb_; pwr(i, j) = pwr_;
                
            end
            
        end
        min(p(:))/2
        findn(p == min(p(:)))
        chi(find(p == min(p(:))))
        
        [ mean(mean(T1_prediction_set, 1)), std(mean(T1_prediction_set, 1)); ...
          mean(mean(FLAIR_prediction_set, 1)), std(mean(FLAIR_prediction_set, 1)); ...
          mean(mean(DTI_prediction_set, 1)), std(mean(DTI_prediction_set, 1)); ...
          mean(mean(fMRI_prediction_set, 1)), std(mean(fMRI_prediction_set, 1)); ...
          mean(mean(mMGlobalf_prediction_set, 1)), std(mean(mMGlobalf_prediction_set, 1)); ...
          mean(mean(mMmS_prediction_set, 1)), std(mean(mMmS_prediction_set, 1)); ...
          mean(mean(mMmSPeri_prediction_set, 1)), std(mean(mMmSPeri_prediction_set, 1)) ]
      
        %% comparison matrix
        %%               T1     FLAIR     DTI     fMRI    mMGlobal  mMGlobalf mMmS mMmSPeri 
        %% T1
        %% FLAIR
        %% DTI
        %% mMGlobal
        %% mMGlobalf
        %% mMmS        0.002   0.038    0.002    0.023     0.007      0.007
        %% mMmSPeri    0.002   0.022  >0.0001    0.013     0.007      0.004
        
    end
        
    %% 7. Do post-evaluation for consistently missed cases (>50%)
    for missing_case = 1
        
        %% 071 (TypeIIA) -> TypeIIB | 076 (TypeIIB) -> control | 078 (TypeIIA) -> control
        table((1:23)', case_num_pat_fMRI, group_class(FCD_idx), [ (mean(mMmSPeri_prediction_set(FCD_idx, :), 2)+0.001).*(mean(mMmSPeri_prediction_set(FCD_idx, :), 2)<0.5) (mode(mMmSPeri_prediction_real_set(FCD_idx, :), 2)).*(mean(mMmSPeri_prediction_set(FCD_idx, :), 2)<0.5) ])
        % [ 11 14 16 ]
        mean_z_lesion_temp(:, [11 14 16]+29)
        
        tabulate(mMmSPeri_prediction_real_set(11+29, :))
        tabulate(mMmSPeri_prediction_real_set(14+29, :))
        tabulate(mMmSPeri_prediction_real_set(16+29, :))
        
        mean_z_lesion_temp(:, [11 14 16]+29)
        
    end  
    
    %% 8. Visualize the classification
    for visualizatoin_using_MDS = 1
          
        mMmSPeri_median_idx = find(abs(mean(mMmSPeri_prediction_set, 1)-0.9231)<0.01);
        discriminative_features = find(mean(squeeze(sum(mMmSPeri_fs_set, 1)>0)') >= 0.3);
        
        %% chose XXth iteration for visualization
        [ [0 mMmSPeri_median_idx]; (1:23)' mMmSPeri_prediction_set(FCD_idx, mMmSPeri_median_idx) ]
        
        i = 57;
        classified_group = group_class;
        classified_group(mMmSPeri_prediction_real_set(:, i)==1) = { 'control' };
        classified_group(mMmSPeri_prediction_real_set(:, i)==2) = { 'FCDIIa' };
        classified_group(mMmSPeri_prediction_real_set(:, i)==3) = { 'FCDIIb' };
        misclassified = (mMmSPeri_prediction_set(:, i)==0);
        % ================================
        % ===== Showing the results ======
        % ================================
        case_num_pat_temp = case_num_pat_fMRI;
        
        % Assign color for each class
        colorList = [ 0 0.8 0; 1 0.5 0.5; 0 0 0 ];
        
                               
        for mds = 1
            
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
            y = newCoor(:,2); y(39) =  0.9;
            patchSize = 30; 
            colorTrueClassPlot = colorTrueClass;
            colorResultClassPlot = colorResultClass;
            figure; S1 = axis; hold on;      
            scatter(x,y,patchSize,colorTrueClassPlot,'filled'); 
            scatter(x(~misclassified),y(~misclassified),patchSize*2,colorResultClass(~misclassified, :)); 
            scatter(x(misclassified),y(misclassified),patchSize*3,[0 0 0], 'x');
            xlim([min(x)-0.5 max(x)+0.5]); ylim([min(x)-0.5 max(x)+0.5]); title('Multidimensional scaling');
            xlabel('x component'); ylabel('y component');  
            
            % legend
            scatter(-3,4,patchSize,colorList(1, :),'filled'); text(-2.5, 4, 'Controls');
            scatter(-3,3.7,patchSize,colorList(2, :),'filled'); text(-2.5, 3.7, 'IIA');
            scatter(-3,3.4,  patchSize,colorList(3, :),'filled'); text(-2.5, 3.4, 'IIB');            
            scatter(-3,3.1,  patchSize*3, [0,0,0], 'x'); text(-2.5, 3.1, 'misclassifed');            
            
            if(save_fig)
                export_fig([OUTPATH '/prediction_multiclass_mds_summary_final_multiclass' ], '-m4', '-png'); close(gcf);
            end
            
        end        
        
    end
    
    %% 9. Permutation test for checking whether it is over the chance level
    for statistical_test_using_mMmSPeri_features = 1
        
        delete(gcp);
        parpool(24);
                
        permutation = 1000;
        iteration_num = 100;
        bagging_num = 10;
               
        opts = [];
        mean_z_lesion_final = mean_z_lesion_temp2;
        keepout_vec = logical(zeros(size(mean_z_lesion_final, 1), 1));
        keepout_vec([ 1 6 10 11 23:28 29 34 38 39 48:53 54 55 57 59 61 63:65 66 67 69 71 73 75:77 80 81 ]) = 1;
        
        prediction_accuracy_rand = zeros(length(group_class), permutation);
        prediction_real_set_rand = zeros(length(group_class), permutation);
        fs_set_rand = zeros(kfold, size(mean_z_lesion_final, 1), permutation);
        
        for p = 1 : permutation
            
            rand_idx        = randperm(length(group_class));
            rand_idx_fold   = randperm(iteration_num); rand_idx_fold = rand_idx_fold(1);  
            histo_type_rand = group_class(rand_idx);
            
            classified_group = cell(length(histo_type_rand), 1);
            
            Indices = kfold_set(rand_idx_fold, :)';
            prediction_curr = zeros(length(histo_type_rand), 1);
            prediction_curr_real = zeros(length(histo_type_rand), 1);
           
            fs_set_temp = cell(kfold, 1);
            history_set_temp = cell(kfold, 1);
            p_val_temp = cell(kfold, 1);
            
            parfor c = 1 : kfold
                
                test_cases = Indices==c;
                CV_group = histo_type_rand;
                CV_group(test_cases) = [];
                mean_z_lesion_temp_lvo = mean_z_lesion_final;
                mean_z_lesion_temp_lvo(:, test_cases) = [];
                
                p = cvpartition(CV_group,'resubstitution');
                [fs,history] = sequentialfs(@fs_func_multiclass_svm, mean_z_lesion_temp_lvo', CV_group,'cv', p, 'direction', 'forward', 'keepout', keepout_vec', 'options',opts);
                fs_set_temp{c} = fs;
                
            end
            fs_set_temp = cell2mat(fs_set_temp);
            
            for c = 1 : kfold
                
                test_cases  = Indices==c;
                train_cases = Indices~=c;
                
                train_control = find(strcmp(histo_type_rand, 'control').*train_cases);
                train_groupA = find(strcmp(histo_type_rand, 'FCDIIa').*train_cases);
                train_groupB = find(strcmp(histo_type_rand, 'FCDIIb').*train_cases);
                
                num_of_control_in_train = length(train_control);
                num_of_FCDIIA_in_train = length(train_groupA);
                num_of_FCDIIB_in_train = length(train_groupB);
                num_of_each_group = min( [ num_of_control_in_train num_of_FCDIIA_in_train num_of_FCDIIB_in_train] );
                
                CV_group = histo_type_rand;
                CV_group(test_cases) = [];
                
                mean_z_lesion_temp_lvo = mean_z_lesion_final;
                mean_z_lesion_temp_lvo(:, test_cases) = [];
                
                test_real  = find(test_cases);
                CV_group_test = histo_type_rand(test_real);
                
                fs = fs_set_temp(c, :);
                
                if(bagging_num==0)
                    
                    classified_group = cell(length(histo_type_rand), 1);
                    train_real = find(train_cases);
                    test_real  = find(test_cases);
                    
                    CV_group_train = histo_type_rand(train_real);
                    CV_group_test = histo_type_rand(test_real);
                    
                    training_data = mean_z_lesion_final(fs, train_real);
                    test_data     = mean_z_lesion_final(fs, test_real);
                    
                    t = templateSVM('KernelFunction', 'polynomial');
                    classOrder = unique(histo_type_rand);
                    
                    Mdl = fitcecoc(training_data',CV_group_train, 'Learners', t, 'ClassNames', classOrder);
                    newClass = predict(Mdl, test_data');
                    
                    ia = [];
                    ia(strcmp(newClass, 'control'), 1) = 1;
                    ia(strcmp(newClass, 'FCDIIa'), 1) = 2;
                    ia(strcmp(newClass, 'FCDIIb'), 1) = 3;
                    
                    prediction_curr(test_cases) = strcmp(CV_group_test, newClass);
                    prediction_curr_real(test_cases) = ia;
                    
                elseif(bagging_num>0)
                    
                    test_bagging = cell(1, bagging_num);
                    CV_group_test = histo_type_rand(test_real);
                    
                    parfor j = 1 : bagging_num
                        
                        train_control_real = sort(train_control(randperm(num_of_control_in_train, num_of_each_group)));
                        train_groupA_real  = sort(train_groupA(randperm(num_of_FCDIIA_in_train, num_of_each_group)));
                        train_groupB_real  = sort(train_groupB(randperm(num_of_FCDIIB_in_train, num_of_each_group)));
                        
                        train_real = sort([ train_control_real; train_groupA_real; train_groupB_real]);
                        test_real  = find(test_cases);
                        
                        CV_group_train = histo_type_rand(train_real);
                        
                        training_data = mean_z_lesion_final(fs, train_real);
                        test_data     = mean_z_lesion_final(fs, test_real);
                        
                        t = templateSVM('KernelFunction', 'polynomial');
                        classOrder = unique(histo_type_rand);
                        
                        Mdl = fitcecoc(training_data',CV_group_train, 'Learners', t, 'ClassNames', classOrder);
                        newClass = predict(Mdl, test_data');
                        
                        ib = [];
                        ib(strcmp(newClass, 'control'), 1) = 1;
                        ib(strcmp(newClass, 'FCDIIa'), 1) = 2;
                        ib(strcmp(newClass, 'FCDIIb'), 1) = 3;
                        
                        test_bagging{j} = ib;
                        
                    end
                    
                    ia = [];
                    ia(strcmp(CV_group_test, 'control'), 1) = 1;
                    ia(strcmp(CV_group_test, 'FCDIIa'), 1) = 2;
                    ia(strcmp(CV_group_test, 'FCDIIb'), 1) = 3;
                    
                    prediction_curr(test_cases) = mode(double(cell2mat(test_bagging)), 2) == ia;
                    prediction_curr_real(test_cases) = mode(double(cell2mat(test_bagging)), 2);
                    
                end
                
            end
            
            disp([ 'iter ', num2str(p), ' : ', num2str(sum(prediction_curr)/length(histo_type_rand))]);
            prediction_accuracy_rand(:, p) = prediction_curr;
            prediction_real_set_rand(:, p) = prediction_curr_real; % 1: control, 2: Type IIA, 3: Type IIB
            fs_set_rand(:, :, p) = fs_set_temp;
            
        end
        
        delete(gcp);
        
    end

end

for t1_flair_dti_fMRI = 1
    
    delete(gcp);
    parpool(24);
    
    %% 1. Make the feature vectors, matching the T1, FLAIR and DTI zscore
    for featureset_build = 1
        
        for control_setup = 1
            
            %% included cases
            case_num_cont = case_num_cont_org;
            case_num_cont_dti = case_num_cont_dti_org;
            case_num_cont_fMRI = case_num_cont_fMRI_org;
            pl_feature_set_T1_FLAIR_cont = pl_feature_set_T1_FLAIR_cont_org;
            pl_feature_set_DTI_cont = pl_feature_set_DTI_cont_org;
            pl_feature_set_fMRI_cont = pl_feature_set_fMRI_cont_org;
            mean_z_lesion_T1_FLAIR_cont = mean_z_lesion_T1_FLAIR_cont_org;
            mean_z_lesion_DTI_cont = mean_z_lesion_DTI_cont_org;
            mean_z_lesion_fMRI_cont = mean_z_lesion_fMRI_cont_org;
            
            for included_case = 1
                
                %% excluded some outliers and data cases newly acquired ...
                
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
                
                %% final combination, all parameters setup according to findings (figure 3)
                % T1:       Cortical thickness, up to 2 mm
                % FLAIR:    RI intracortical, up to 6 mm
                % DTI:      Cortical FA, up to 4 mm
                % fMRI:     ReHo, up to 16 mm
                mean_z_perilesion_T1_cont_final    = mean(pl_feature_set_T1_FLAIR_cont(6).data, 4)'; mean_z_perilesion_T1_cont_final = mean(mean_z_perilesion_T1_cont_final(:, 1:3), 2)';
                mean_z_perilesion_FLAIR_cont_final = mean(squeeze(mean(pl_feature_set_T1_FLAIR_cont(9).data(:, :, 2:5, :), 3)), 3)'; mean_z_perilesion_FLAIR_cont_final = mean(mean_z_perilesion_FLAIR_cont_final(:, 1:6), 2)';
                mean_z_perilesion_DTI_cont_final   = mean(squeeze(mean(pl_feature_set_DTI_cont(1).data(:, :, 2:5, :), 3)), 3)'; mean_z_perilesion_DTI_cont_final = mean(mean_z_perilesion_DTI_cont_final(:, 1:5), 2)';
                mean_z_perilesion_fMRI_cont_final  = mean(pl_feature_set_fMRI_cont(2).data, 4)'; mean_z_perilesion_fMRI_cont_final = mean(mean_z_perilesion_fMRI_cont_final(:, 1:14), 2)';
                
                % original combining all lesional and peri-lesional features & one without peri-lesional abnormalities
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
            
            %% excluded cases: 319, 320, 321, 323, 324, 325, 326,327, 329, those who don't have fMRI
            excluded_cases = { '306_1', '313_1', '322_1' };
            [C, ia1, ib1] = intersect(case_num_cont_org, excluded_cases);
            case_num_cont_org_ = case_num_cont_org;
            case_num_cont_org_(ia1) = [];
            [C, ia1, ib1] = intersect(case_num_cont, case_num_cont_org_);
            case_num_cont_exc = case_num_cont_org_;
            case_num_cont_exc(ib1) = [];
            [C, iaexc, ibexc] = intersect(case_num_cont_exc, case_num_cont_org_);
            mean_z_lesion_T1_FLAIR_cont_exc         = mean_z_lesion_T1_FLAIR_cont_org(:, :, ibexc);
            mean_z_lesion_DTI_cont_exc              = mean_z_lesion_DTI_cont_org(:, :, ibexc);
            mean_z_lesion_fMRI_cont_exc             = repmat(mean(mean_z_lesion_fMRI_cont, 3), 1, 1, length(ibexc));
            mean_z_perilesion_T1_cont_exc           = mean(pl_feature_set_T1_FLAIR_cont_org(6).data(peri_lesional_dist, ibexc, :, :), 4);
            mean_z_perilesion_FLAIR_cont_exc        = mean(squeeze(mean(pl_feature_set_T1_FLAIR_cont_org(9).data(peri_lesional_dist, ibexc, 2:5, :), 3)), 3);
            mean_z_perilesion_DTI_cont_exc          = repmat(mean(mean(squeeze(mean(pl_feature_set_DTI_cont(1).data(peri_lesional_dist, :, 2:5, :), 3)), 3), 2), 1, length(ibexc));
            mean_z_perilesion_fMRI_cont_exc1        = repmat(mean(mean(pl_feature_set_fMRI_cont(1).data(peri_lesional_dist_fMRI, :, :, :), 4), 2), 1, length(ibexc));
            mean_z_perilesion_fMRI_cont_exc2        = repmat(mean(mean(pl_feature_set_fMRI_cont(2).data(peri_lesional_dist_fMRI, :, :, :), 4), 2), 1, length(ibexc));
            mean_z_lesion_cont_exc                  = [ squeeze(mean(mean_z_lesion_T1_FLAIR_cont_exc, 2)); squeeze(mean(mean_z_lesion_DTI_cont_exc, 2)); squeeze(mean(mean_z_lesion_fMRI_cont_exc, 2)); ...
                                                        mean_z_perilesion_T1_cont_exc; mean_z_perilesion_FLAIR_cont_exc; mean_z_perilesion_DTI_cont_exc; mean_z_perilesion_fMRI_cont_exc1; mean_z_perilesion_fMRI_cont_exc2 ];
            
        end

        for patient_setup = 1
            
            %% included cases
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
            
            %% original one with periliesional anomalies & one without
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
            
            %% excluded cases
            [C, ia1, ib1] = intersect(case_num_pat_fMRI, case_num_pat);
            case_num_pat_exc = case_num_pat;
            case_num_pat_exc(ib1) = [];
            [C, iaexc, ibexc] = intersect(case_num_pat_exc, case_num_pat);
            histo_type_exc = histo_type(ibexc);
            mean_z_lesion_T1_FLAIR_pat_exc = mean_z_lesion_T1_FLAIR(:, ibexc);
            mean_z_lesion_DTI_pat_exc = zeros(24, length(ibexc));
            mean_z_lesion_DTI_pat_exc(:, strcmp(histo_type(ibexc), 'FCDIIb'))  = repmat(mean(mean_z_lesion_DTI(:, strcmp(histo_type_dti, 'FCDIIb')), 2), 1, sum(strcmp(histo_type(ibexc), 'FCDIIb')));
            mean_z_lesion_DTI_pat_exc(:, strcmp(histo_type(ibexc), 'FCDIIa'))  = repmat(mean(mean_z_lesion_DTI(:, strcmp(histo_type_dti, 'FCDIIa')), 2), 1, sum(strcmp(histo_type(ibexc), 'FCDIIa')));
            mean_z_lesion_fMRI_pat_exc = zeros(4, length(ibexc));
            mean_z_lesion_fMRI_pat_exc(:, strcmp(histo_type(ibexc), 'FCDIIb')) = repmat(mean(mean_z_lesion_fMRI(:, strcmp(histo_type_fMRI, 'FCDIIb')), 2), 1, sum(strcmp(histo_type(ibexc), 'FCDIIb')));
            mean_z_lesion_fMRI_pat_exc(:, strcmp(histo_type(ibexc), 'FCDIIa')) = repmat(mean(mean_z_lesion_fMRI(:, strcmp(histo_type_fMRI, 'FCDIIa')), 2), 1, sum(strcmp(histo_type(ibexc), 'FCDIIa')));
            mean_z_perilesion_T1_exc                                           = pl_feature_set_T1_FLAIR(6).data(peri_lesional_dist, ibexc);
            mean_z_perilesion_FLAIR_exc                                        = mean(pl_feature_set_T1_FLAIR(9).data(peri_lesional_dist, ibexc, 2:5), 3);
            mean_z_perilesion_DTI_exc = zeros(2, length(ibexc));
            mean_z_perilesion_DTI_exc(:, strcmp(histo_type(ibexc), 'FCDIIb'))  = repmat(mean(mean(pl_feature_set_DTI(1).data(peri_lesional_dist, strcmp(histo_type_dti, 'FCDIIb'), [3 5]), 3), 2), 1, sum(strcmp(histo_type(ibexc), 'FCDIIb')));
            mean_z_perilesion_DTI_exc(:, strcmp(histo_type(ibexc), 'FCDIIa'))  = repmat(mean(mean(pl_feature_set_DTI(1).data(peri_lesional_dist, strcmp(histo_type_dti, 'FCDIIa'), [3 5]), 3), 2), 1, sum(strcmp(histo_type(ibexc), 'FCDIIa')));
            mean_z_perilesion_fMRI_exc1 = zeros(2, length(ibexc));
            mean_z_perilesion_fMRI_exc1(:, strcmp(histo_type(ibexc), 'FCDIIb'))  = repmat(mean(pl_feature_set_fMRI(1).data(peri_lesional_dist_fMRI, strcmp(histo_type_fMRI, 'FCDIIb')), 2), 1, sum(strcmp(histo_type(ibexc), 'FCDIIb')));
            mean_z_perilesion_fMRI_exc1(:, strcmp(histo_type(ibexc), 'FCDIIa'))  = repmat(mean(pl_feature_set_fMRI(1).data(peri_lesional_dist_fMRI, strcmp(histo_type_fMRI, 'FCDIIa')), 2), 1, sum(strcmp(histo_type(ibexc), 'FCDIIa')));
            mean_z_perilesion_fMRI_exc2 = zeros(2, length(ibexc));
            mean_z_perilesion_fMRI_exc2(:, strcmp(histo_type(ibexc), 'FCDIIb'))  = repmat(mean(pl_feature_set_fMRI(2).data(peri_lesional_dist_fMRI, strcmp(histo_type_fMRI, 'FCDIIb')), 2), 1, sum(strcmp(histo_type(ibexc), 'FCDIIb')));
            mean_z_perilesion_fMRI_exc2(:, strcmp(histo_type(ibexc), 'FCDIIa'))  = repmat(mean(pl_feature_set_fMRI(2).data(peri_lesional_dist_fMRI, strcmp(histo_type_fMRI, 'FCDIIa')), 2), 1, sum(strcmp(histo_type(ibexc), 'FCDIIa')));
            mean_z_lesion_pat_exc                  = [ mean_z_lesion_T1_FLAIR_pat_exc; mean_z_lesion_DTI_pat_exc; mean_z_lesion_fMRI_pat_exc; ...
                                                       mean_z_perilesion_T1_exc; mean_z_perilesion_FLAIR_exc; mean_z_perilesion_DTI_exc; mean_z_perilesion_fMRI_exc1; mean_z_perilesion_fMRI_exc2 ];                

        end
        
        mean_z_lesion        = [ mean_z_lesion_cont mean_z_lesion_FCD             ];                    %% 85x52: 81(lesional)+4(peri) data with perilesional anomalies, peri-lesioanl being all averaged
        mean_z_lesion_temp   = [ mean_z_lesion_cont_temp mean_z_lesion_FCD_temp   ];                    %% 81x52: 81(lesional) data without perilesional anomalies
        mean_z_lesion_temp2  = [ mean_z_lesion_cont_temp2 mean_z_lesion_FCD_temp2 ];                    %% 91x52: 81(lesional)+10(peri) data with perilesional anomalies, peri-lesioanl being all seperately across distance
        mean_z_lesion_exc    = [ mean_z_lesion_cont_exc mean_z_lesion_pat_exc     ];                    %% 91x19: 81(lesional)+10(peri) x 12(control)+7(FCD) data with perilesional anomalies, peri-lesioanl being all seperately across distance
        
        group_class = cell(size(mean_z_lesion_cont, 2)+size(mean_z_lesion_FCD, 2), 1);
        group_class(1:size(mean_z_lesion_cont, 2)) = { 'control' };
        group_class(size(mean_z_lesion_cont, 2)+1:end) = histo_type_fMRI;
        
        group_class_exc = cell(size(mean_z_lesion_cont_exc, 2)+size(mean_z_lesion_pat_exc, 2), 1);
        group_class_exc(1:size(mean_z_lesion_cont_exc, 2)) = { 'control' };
        group_class_exc(size(mean_z_lesion_cont_exc, 2)+1:end) = histo_type_exc;
        
    end
    
    %% 2. Predefine k-fold cross-validation sets
    iteration_num = 100;     
    
    for i = 1 : iteration_num
        
        Indices = crossvalind('kfold', group_class, kfold);
        kfold_set(i, :) = Indices';        
        
    end
    
    bagging_num = 100;

    %% 3. Use unimodal global features based on pre-selection
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
        %% 82 83 = Thickness of T1 perilesion    -> 82 83    p24
        %% 84 85 = Intensity of FLAIR perilesion -> 84 85    p25
        %% 86 87 = FA of DTI perilesion          -> 86 87    p26
        %% 88 89 = ALFF perilesion               -> 88 89    p27
        %% 90 91 = ReHo perilesion               -> 90 91    p28
        
        %% Using T1 global
        for T1_glboal_features = 1
            
            % one vs. one 
            warning off;
            classified_group = group_class;
            SVMstruct_set = cell(length(classified_group), 1);
            discriminative_features_set = cell(length(group_class), 1);      
            discriminative_features = { [2:5], [7:9 16], [12:15], 17, 18, 19, [20:22] };
            opts = statset();      
            
            mean_z_lesion_final = [];
            for i = 1 : length(discriminative_features)
                mean_z_lesion_final = [ mean_z_lesion_final; mean(mean_z_lesion_temp(discriminative_features{i}, :), 1) ];
            end
                        
            keepout_vec = logical(zeros(size(mean_z_lesion_final, 1), 1));
            keepout_vec([  ]) = 1;
            
            CONTROL    = strcmp(group_class, 'control');
            FCDTypeIIA = strcmp(group_class, 'FCDIIa');
            FCDTypeIIB = strcmp(group_class, 'FCDIIb');
            
            prediction_set = zeros(length(group_class), iteration_num);
            prediction_real_set = zeros(length(group_class), iteration_num);
            fs_set = zeros(kfold, size(mean_z_lesion_final, 1), iteration_num);
            
            for i = 1 : iteration_num
                
                Indices = kfold_set(i, :)';
                prediction_curr = zeros(length(group_class), 1);
                prediction_curr_real = zeros(length(group_class), 1);
                
                fs_set_temp = cell(kfold, 1);
                history_set_temp = cell(kfold, 1);
                p_val_temp = cell(kfold, 1);
                
                parfor c = 1 : kfold
                    
                    test_cases = Indices==c;                    
                    
                    CV_group = group_class;
                    CV_group(test_cases) = [];
                    
                    mean_z_lesion_temp_lvo = mean_z_lesion_final;
                    mean_z_lesion_temp_lvo(:, test_cases) = [];
                    
                    p = cvpartition(CV_group,'resubstitution');
                    [fs,history] = sequentialfs(@fs_func_multiclass_svm, mean_z_lesion_temp_lvo', CV_group,'cv', p, 'direction', 'forward', 'keepout', keepout_vec', 'options',opts);
                    fs_set_temp{c} = fs;
                    
                end
                
                fs_set_temp = cell2mat(fs_set_temp);                
                
                for c = 1 : kfold
                    
                    test_cases  = Indices==c;
                    train_cases = Indices~=c;
                    
                    train_control = find(strcmp(group_class, 'control').*train_cases);
                    train_groupA = find(strcmp(group_class, 'FCDIIa').*train_cases);
                    train_groupB = find(strcmp(group_class, 'FCDIIb').*train_cases);
                    
                    num_of_control_in_train = length(train_control);
                    num_of_FCDIIA_in_train = length(train_groupA);
                    num_of_FCDIIB_in_train = length(train_groupB);
                    num_of_each_group = min( [ num_of_control_in_train num_of_FCDIIA_in_train num_of_FCDIIB_in_train] );
                    
                    CV_group = group_class;
                    CV_group(test_cases) = [];
                    
                    mean_z_lesion_temp_lvo = mean_z_lesion_final;
                    mean_z_lesion_temp_lvo(:, test_cases) = [];
                    
                    test_real  = find(test_cases);
                    CV_group_test = group_class(test_real);
                    
                    fs = fs_set_temp(c, :);
                    
                    if(bagging_num==0)
                        
                        classified_group = cell(length(group_class), 1);
                        train_real = find(train_cases);
                        test_real  = find(test_cases);
                        
                        CV_group_train = group_class(train_real);
                        CV_group_test = group_class(test_real);
                        
                        training_data = mean_z_lesion_final(fs, train_real);
                        test_data     = mean_z_lesion_final(fs, test_real);
                        
                        t = templateSVM('KernelFunction', 'polynomial');
                        classOrder = unique(group_class);
                        
                        Mdl = fitcecoc(training_data',CV_group_train, 'Learners', t, 'ClassNames', classOrder);
                        newClass = predict(Mdl, test_data');
                        
                        ia = [];
                        ia(strcmp(newClass, 'control'), 1) = 1;
                        ia(strcmp(newClass, 'FCDIIa'), 1) = 2;
                        ia(strcmp(newClass, 'FCDIIb'), 1) = 3;
                        
                        prediction_curr(test_cases) = strcmp(CV_group_test, newClass);
                        prediction_curr_real(test_cases) = ia;
                        
                    elseif(bagging_num>0)
                        
                        test_bagging = cell(1, bagging_num);
                        CV_group_test = group_class(test_real);
                        
                        parfor j = 1 : bagging_num
                            
                            train_control_real = sort(train_control(randperm(num_of_control_in_train, num_of_each_group)));
                            train_groupA_real  = sort(train_groupA(randperm(num_of_FCDIIA_in_train, num_of_each_group)));
                            train_groupB_real  = sort(train_groupB(randperm(num_of_FCDIIB_in_train, num_of_each_group)));
                            
                            train_real = sort([ train_control_real; train_groupA_real; train_groupB_real]);
                            test_real  = find(test_cases);
                            
                            CV_group_train = group_class(train_real);                            
                            
                            training_data = mean_z_lesion_final(fs, train_real);
                            test_data     = mean_z_lesion_final(fs, test_real);
                            
                            t = templateSVM('KernelFunction', 'polynomial');
                            classOrder = unique(group_class);
                            
                            Mdl = fitcecoc(training_data',CV_group_train, 'Learners', t, 'ClassNames', classOrder);
                            newClass = predict(Mdl, test_data');                            
                            
                            ib = [];
                            ib(strcmp(newClass, 'control'), 1) = 1;
                            ib(strcmp(newClass, 'FCDIIa'), 1) = 2;
                            ib(strcmp(newClass, 'FCDIIb'), 1) = 3;
                            
                            test_bagging{j} = ib;
                            
                        end
                        
                        ia = [];
                        ia(strcmp(CV_group_test, 'control'), 1) = 1;
                        ia(strcmp(CV_group_test, 'FCDIIa'), 1) = 2;
                        ia(strcmp(CV_group_test, 'FCDIIb'), 1) = 3;
                        
                        prediction_curr(test_cases) = mode(double(cell2mat(test_bagging)), 2) == ia;
                        prediction_curr_real(test_cases) = mode(double(cell2mat(test_bagging)), 2);
                        
                    end
                    
                end
                    
                disp([ 'iter ', num2str(i), ' : ', num2str(sum(prediction_curr)/length(group_class))]);
                prediction_set(:, i) = prediction_curr;
                prediction_real_set(:, i) = prediction_curr_real; % 1: control, 2: Type IIA, 3: Type IIB
                fs_set(:, :, i) = fs_set_temp;
                
            end
            
            T1_prediction_set = prediction_set;
            T1_prediction_real_set = prediction_real_set;
            T1_fs_set = fs_set;
            
        end
        
        %% Using FLAIR global
        for FLAIR_glboal_features = 1
                        
            % one vs. one 
            warning off;
            classified_group = group_class;
            SVMstruct_set = cell(length(classified_group), 1);
            discriminative_features_set = cell(length(group_class), 1);      
            discriminative_features = { [30:33], [35:37 44], [40:43], [45:47] };
            opts = statset();      
            
            mean_z_lesion_final = [];
            for i = 1 : length(discriminative_features)
                mean_z_lesion_final = [ mean_z_lesion_final; mean(mean_z_lesion_temp(discriminative_features{i}, :), 1) ];
            end
                        
            keepout_vec = logical(zeros(size(mean_z_lesion_final, 1), 1));
            keepout_vec([ ]) = 1;
            
            CONTROL    = strcmp(group_class, 'control');
            FCDTypeIIA = strcmp(group_class, 'FCDIIa');
            FCDTypeIIB = strcmp(group_class, 'FCDIIb');
            
            prediction_set = zeros(length(group_class), iteration_num);
            prediction_real_set = zeros(length(group_class), iteration_num);
            fs_set = zeros(kfold, size(mean_z_lesion_final, 1), iteration_num);            
            
            for i = 1 : iteration_num
                
                Indices = kfold_set(i, :)';
                prediction_curr = zeros(length(group_class), 1);
                prediction_curr_real = zeros(length(group_class), 1);
                
                fs_set_temp = cell(kfold, 1);
                history_set_temp = cell(kfold, 1);
                p_val_temp = cell(kfold, 1);
                
                parfor c = 1 : kfold
                    
                    test_cases = Indices==c;                    
                    
                    CV_group = group_class;
                    CV_group(test_cases) = [];
                    
                    mean_z_lesion_temp_lvo = mean_z_lesion_final;
                    mean_z_lesion_temp_lvo(:, test_cases) = [];
                    
                    p = cvpartition(CV_group,'resubstitution');
                    [fs,history] = sequentialfs(@fs_func_multiclass_svm, mean_z_lesion_temp_lvo', CV_group,'cv', p, 'direction', 'forward', 'keepout', keepout_vec', 'options',opts);
                    fs_set_temp{c} = fs;
                    
                end
                
                fs_set_temp = cell2mat(fs_set_temp);                
                
                for c = 1 : kfold
                    
                    test_cases  = Indices==c;
                    train_cases = Indices~=c;
                    
                    train_control = find(strcmp(group_class, 'control').*train_cases);
                    train_groupA = find(strcmp(group_class, 'FCDIIa').*train_cases);
                    train_groupB = find(strcmp(group_class, 'FCDIIb').*train_cases);
                    
                    num_of_control_in_train = length(train_control);
                    num_of_FCDIIA_in_train = length(train_groupA);
                    num_of_FCDIIB_in_train = length(train_groupB);
                    num_of_each_group = min( [ num_of_control_in_train num_of_FCDIIA_in_train num_of_FCDIIB_in_train] );
                    
                    CV_group = group_class;
                    CV_group(test_cases) = [];
                    
                    mean_z_lesion_temp_lvo = mean_z_lesion_final;
                    mean_z_lesion_temp_lvo(:, test_cases) = [];
                    
                    test_real  = find(test_cases);
                    CV_group_test = group_class(test_real);
                    
                    fs = fs_set_temp(c, :);
                    
                    if(bagging_num==0)
                        
                        classified_group = cell(length(group_class), 1);
                        train_real = find(train_cases);
                        test_real  = find(test_cases);
                        
                        CV_group_train = group_class(train_real);
                        CV_group_test = group_class(test_real);
                        
                        training_data = mean_z_lesion_final(fs, train_real);
                        test_data     = mean_z_lesion_final(fs, test_real);
                        
                        t = templateSVM('KernelFunction', 'polynomial', 'Standardize', 1, 'KernelScale', 'auto');
                        classOrder = unique(group_class);
                        
                        Mdl = fitcecoc(training_data',CV_group_train, 'Learners', t, 'ClassNames', classOrder);
                        newClass = predict(Mdl, test_data');
                        
                        ia = [];
                        ia(strcmp(newClass, 'control'), 1) = 1;
                        ia(strcmp(newClass, 'FCDIIa'), 1) = 2;
                        ia(strcmp(newClass, 'FCDIIb'), 1) = 3;
                        
                        prediction_curr(test_cases) = strcmp(CV_group_test, newClass);
                        prediction_curr_real(test_cases) = ia;
                        
                    elseif(bagging_num>0)
                        
                        test_bagging = cell(1, bagging_num);
                        CV_group_test = group_class(test_real);
                        
                        parfor j = 1 : bagging_num
                            
                            train_control_real = sort(train_control(randperm(num_of_control_in_train, num_of_each_group)));
                            train_groupA_real  = sort(train_groupA(randperm(num_of_FCDIIA_in_train, num_of_each_group)));
                            train_groupB_real  = sort(train_groupB(randperm(num_of_FCDIIB_in_train, num_of_each_group)));
                            
                            train_real = sort([ train_control_real; train_groupA_real; train_groupB_real]);
                            test_real  = find(test_cases);
                            
                            CV_group_train = group_class(train_real);                            
                            
                            training_data = mean_z_lesion_final(fs, train_real);
                            test_data     = mean_z_lesion_final(fs, test_real);
                            
                            t = templateSVM('KernelFunction', 'polynomial', 'Standardize', 1, 'KernelScale', 'auto');
                            classOrder = unique(group_class);
                            
                            Mdl = fitcecoc(training_data',CV_group_train, 'Learners', t, 'ClassNames', classOrder);
                            newClass = predict(Mdl, test_data');                            
                            
                            ib = [];
                            ib(strcmp(newClass, 'control'), 1) = 1;
                            ib(strcmp(newClass, 'FCDIIa'), 1) = 2;
                            ib(strcmp(newClass, 'FCDIIb'), 1) = 3;
                            
                            test_bagging{j} = ib;
                            
                        end
                        
                        ia = [];
                        ia(strcmp(CV_group_test, 'control'), 1) = 1;
                        ia(strcmp(CV_group_test, 'FCDIIa'), 1) = 2;
                        ia(strcmp(CV_group_test, 'FCDIIb'), 1) = 3;
                        
                        prediction_curr(test_cases) = mode(double(cell2mat(test_bagging)), 2) == ia;
                        prediction_curr_real(test_cases) = mode(double(cell2mat(test_bagging)), 2);
                        
                    end
                    
                end
                    
                disp([ 'iter ', num2str(i), ' : ', num2str(sum(prediction_curr)/length(group_class))]);
                prediction_set(:, i) = prediction_curr;
                prediction_real_set(:, i) = prediction_curr_real; % 1: control, 2: Type IIA, 3: Type IIB
                fs_set(:, :, i) = fs_set_temp;
                
            end
            
            FLAIR_prediction_set = prediction_set;
            FLAIR_prediction_real_set = prediction_real_set;
            FLAIR_fs_set = fs_set;
            
        end
        
        %% Using DTI global
        for DTI_glboal_features = 1
                      
            % one vs. one 
            warning off;
            classified_group = group_class;
            SVMstruct_set = cell(length(classified_group), 1);
            discriminative_features_set = cell(length(group_class), 1);      
            discriminative_features = { [56 58], [60 62], [68 70], [72 74] };
            opts = statset();      
            
            mean_z_lesion_final = [];
            for i = 1 : length(discriminative_features)
                mean_z_lesion_final = [ mean_z_lesion_final; mean(mean_z_lesion_temp(discriminative_features{i}, :), 1) ];
            end
                        
            keepout_vec = logical(zeros(size(mean_z_lesion_final, 1), 1));
            keepout_vec([  ]) = 1;
            
            CONTROL    = strcmp(group_class, 'control');
            FCDTypeIIA = strcmp(group_class, 'FCDIIa');
            FCDTypeIIB = strcmp(group_class, 'FCDIIb');
            
            prediction_set = zeros(length(group_class), iteration_num);
            prediction_real_set = zeros(length(group_class), iteration_num);
            fs_set = zeros(kfold, size(mean_z_lesion_final, 1), iteration_num);
            
            for i = 1 : iteration_num
                
                Indices = kfold_set(i, :)';
                prediction_curr = zeros(length(group_class), 1);
                prediction_curr_real = zeros(length(group_class), 1);
                
                fs_set_temp = cell(kfold, 1);
                history_set_temp = cell(kfold, 1);
                p_val_temp = cell(kfold, 1);
                
                parfor c = 1 : kfold
                    
                    test_cases = Indices==c;                    
                    
                    CV_group = group_class;
                    CV_group(test_cases) = [];
                    
                    mean_z_lesion_temp_lvo = mean_z_lesion_final;
                    mean_z_lesion_temp_lvo(:, test_cases) = [];
                    
                    p = cvpartition(CV_group,'resubstitution');
                    [fs,history] = sequentialfs(@fs_func_multiclass_svm, mean_z_lesion_temp_lvo', CV_group,'cv', p, 'direction', 'forward', 'keepout', keepout_vec', 'options',opts);
                    fs_set_temp{c} = fs;
                    
                end
                
                fs_set_temp = cell2mat(fs_set_temp);                
                
                for c = 1 : kfold
                    
                    test_cases  = Indices==c;
                    train_cases = Indices~=c;
                    
                    train_control = find(strcmp(group_class, 'control').*train_cases);
                    train_groupA = find(strcmp(group_class, 'FCDIIa').*train_cases);
                    train_groupB = find(strcmp(group_class, 'FCDIIb').*train_cases);
                    
                    num_of_control_in_train = length(train_control);
                    num_of_FCDIIA_in_train = length(train_groupA);
                    num_of_FCDIIB_in_train = length(train_groupB);
                    num_of_each_group = min( [ num_of_control_in_train num_of_FCDIIA_in_train num_of_FCDIIB_in_train] );
                    
                    CV_group = group_class;
                    CV_group(test_cases) = [];
                    
                    mean_z_lesion_temp_lvo = mean_z_lesion_final;
                    mean_z_lesion_temp_lvo(:, test_cases) = [];
                    
                    test_real  = find(test_cases);
                    CV_group_test = group_class(test_real);
                    
                    fs = fs_set_temp(c, :);
                    
                    if(bagging_num==0)
                        
                        classified_group = cell(length(group_class), 1);
                        train_real = find(train_cases);
                        test_real  = find(test_cases);
                        
                        CV_group_train = group_class(train_real);
                        CV_group_test = group_class(test_real);
                        
                        training_data = mean_z_lesion_final(fs, train_real);
                        test_data     = mean_z_lesion_final(fs, test_real);
                        
                        t = templateSVM('KernelFunction', 'polynomial');
                        classOrder = unique(group_class);
                        
                        Mdl = fitcecoc(training_data',CV_group_train, 'Learners', t, 'ClassNames', classOrder);
                        newClass = predict(Mdl, test_data');
                        
                        ia = [];
                        ia(strcmp(newClass, 'control'), 1) = 1;
                        ia(strcmp(newClass, 'FCDIIa'), 1) = 2;
                        ia(strcmp(newClass, 'FCDIIb'), 1) = 3;
                        
                        prediction_curr(test_cases) = strcmp(CV_group_test, newClass);
                        prediction_curr_real(test_cases) = ia;
                        
                    elseif(bagging_num>0)
                        
                        test_bagging = cell(1, bagging_num);
                        CV_group_test = group_class(test_real);
                        
                        parfor j = 1 : bagging_num
                            
                            train_control_real = sort(train_control(randperm(num_of_control_in_train, num_of_each_group)));
                            train_groupA_real  = sort(train_groupA(randperm(num_of_FCDIIA_in_train, num_of_each_group)));
                            train_groupB_real  = sort(train_groupB(randperm(num_of_FCDIIB_in_train, num_of_each_group)));
                            
                            train_real = sort([ train_control_real; train_groupA_real; train_groupB_real]);
                            test_real  = find(test_cases);
                            
                            CV_group_train = group_class(train_real);                            
                            
                            training_data = mean_z_lesion_final(fs, train_real);
                            test_data     = mean_z_lesion_final(fs, test_real);
                            
                            t = templateSVM('KernelFunction', 'polynomial');
                            classOrder = unique(group_class);
                            
                            Mdl = fitcecoc(training_data',CV_group_train, 'Learners', t, 'ClassNames', classOrder);
                            newClass = predict(Mdl, test_data');                            
                            
                            ib = [];
                            ib(strcmp(newClass, 'control'), 1) = 1;
                            ib(strcmp(newClass, 'FCDIIa'), 1) = 2;
                            ib(strcmp(newClass, 'FCDIIb'), 1) = 3;
                            
                            test_bagging{j} = ib;
                            
                        end
                        
                        ia = [];
                        ia(strcmp(CV_group_test, 'control'), 1) = 1;
                        ia(strcmp(CV_group_test, 'FCDIIa'), 1) = 2;
                        ia(strcmp(CV_group_test, 'FCDIIb'), 1) = 3;
                        
                        prediction_curr(test_cases) = mode(double(cell2mat(test_bagging)), 2) == ia;
                        prediction_curr_real(test_cases) = mode(double(cell2mat(test_bagging)), 2);
                        
                    end
                    
                end
                    
                disp([ 'iter ', num2str(i), ' : ', num2str(sum(prediction_curr)/length(group_class))]);
                prediction_set(:, i) = prediction_curr;
                prediction_real_set(:, i) = prediction_curr_real; % 1: control, 2: Type IIA, 3: Type IIB
                fs_set(:, :, i) = fs_set_temp;
                
            end
            
            DTI_prediction_set = prediction_set;
            DTI_prediction_real_set = prediction_real_set;
            DTI_fs_set = fs_set;
            
        end
        
        %% Using functional global
        for functional_glboal_features = 1
                        
            % one vs. one 
            warning off;
            classified_group = group_class;
            SVMstruct_set = cell(length(classified_group), 1);
            discriminative_features_set = cell(length(group_class), 1);      
            discriminative_features = { 78, 79 };
            opts = statset();      
            
            mean_z_lesion_final = [];
            for i = 1 : length(discriminative_features)
                mean_z_lesion_final = [ mean_z_lesion_final; mean(mean_z_lesion_temp(discriminative_features{i}, :), 1) ];
            end
                        
            keepout_vec = logical(zeros(size(mean_z_lesion_final, 1), 1));
            keepout_vec([  ]) = 1;
            
            CONTROL    = strcmp(group_class, 'control');
            FCDTypeIIA = strcmp(group_class, 'FCDIIa');
            FCDTypeIIB = strcmp(group_class, 'FCDIIb');
            
            prediction_set = zeros(length(group_class), iteration_num);
            prediction_real_set = zeros(length(group_class), iteration_num);
            fs_set = zeros(kfold, size(mean_z_lesion_final, 1), iteration_num);
            
            for i = 1 : iteration_num
                
                Indices = kfold_set(i, :)';
                prediction_curr = zeros(length(group_class), 1);
                prediction_curr_real = zeros(length(group_class), 1);
                
                fs_set_temp = cell(kfold, 1);
                history_set_temp = cell(kfold, 1);
                p_val_temp = cell(kfold, 1);
                
                parfor c = 1 : kfold
                    
                    test_cases = Indices==c;                    
                    
                    CV_group = group_class;
                    CV_group(test_cases) = [];
                    
                    mean_z_lesion_temp_lvo = mean_z_lesion_final;
                    mean_z_lesion_temp_lvo(:, test_cases) = [];
                    
                    p = cvpartition(CV_group,'resubstitution');
                    [fs,history] = sequentialfs(@fs_func_multiclass_svm, mean_z_lesion_temp_lvo', CV_group,'cv', p, 'direction', 'forward', 'keepout', keepout_vec', 'options',opts);
                    fs_set_temp{c} = fs;
                    
                end
                
                fs_set_temp = cell2mat(fs_set_temp);                
                
                for c = 1 : kfold
                    
                    test_cases  = Indices==c;
                    train_cases = Indices~=c;
                    
                    train_control = find(strcmp(group_class, 'control').*train_cases);
                    train_groupA = find(strcmp(group_class, 'FCDIIa').*train_cases);
                    train_groupB = find(strcmp(group_class, 'FCDIIb').*train_cases);
                    
                    num_of_control_in_train = length(train_control);
                    num_of_FCDIIA_in_train = length(train_groupA);
                    num_of_FCDIIB_in_train = length(train_groupB);
                    num_of_each_group = min( [ num_of_control_in_train num_of_FCDIIA_in_train num_of_FCDIIB_in_train] );
                    
                    CV_group = group_class;
                    CV_group(test_cases) = [];
                    
                    mean_z_lesion_temp_lvo = mean_z_lesion_final;
                    mean_z_lesion_temp_lvo(:, test_cases) = [];
                    
                    test_real  = find(test_cases);
                    CV_group_test = group_class(test_real);
                    
                    fs = fs_set_temp(c, :);
                    
                    if(bagging_num==0)
                        
                        classified_group = cell(length(group_class), 1);
                        train_real = find(train_cases);
                        test_real  = find(test_cases);
                        
                        CV_group_train = group_class(train_real);
                        CV_group_test = group_class(test_real);
                        
                        training_data = mean_z_lesion_final(fs, train_real);
                        test_data     = mean_z_lesion_final(fs, test_real);
                        
                        t = templateSVM('KernelFunction', 'polynomial');
                        classOrder = unique(group_class);
                        
                        Mdl = fitcecoc(training_data',CV_group_train, 'Learners', t, 'ClassNames', classOrder);
                        newClass = predict(Mdl, test_data');
                        
                        ia = [];
                        ia(strcmp(newClass, 'control'), 1) = 1;
                        ia(strcmp(newClass, 'FCDIIa'), 1) = 2;
                        ia(strcmp(newClass, 'FCDIIb'), 1) = 3;
                        
                        prediction_curr(test_cases) = strcmp(CV_group_test, newClass);
                        prediction_curr_real(test_cases) = ia;
                        
                    elseif(bagging_num>0)
                        
                        test_bagging = cell(1, bagging_num);
                        CV_group_test = group_class(test_real);
                        
                        parfor j = 1 : bagging_num
                            
                            train_control_real = sort(train_control(randperm(num_of_control_in_train, num_of_each_group)));
                            train_groupA_real  = sort(train_groupA(randperm(num_of_FCDIIA_in_train, num_of_each_group)));
                            train_groupB_real  = sort(train_groupB(randperm(num_of_FCDIIB_in_train, num_of_each_group)));
                            
                            train_real = sort([ train_control_real; train_groupA_real; train_groupB_real]);
                            test_real  = find(test_cases);
                            
                            CV_group_train = group_class(train_real);                            
                            
                            training_data = mean_z_lesion_final(fs, train_real);
                            test_data     = mean_z_lesion_final(fs, test_real);
                            
                            t = templateSVM('KernelFunction', 'polynomial');
                            classOrder = unique(group_class);
                            
                            Mdl = fitcecoc(training_data',CV_group_train, 'Learners', t, 'ClassNames', classOrder);
                            newClass = predict(Mdl, test_data');                            
                            
                            ib = [];
                            ib(strcmp(newClass, 'control'), 1) = 1;
                            ib(strcmp(newClass, 'FCDIIa'), 1) = 2;
                            ib(strcmp(newClass, 'FCDIIb'), 1) = 3;
                            
                            test_bagging{j} = ib;
                            
                        end
                        
                        ia = [];
                        ia(strcmp(CV_group_test, 'control'), 1) = 1;
                        ia(strcmp(CV_group_test, 'FCDIIa'), 1) = 2;
                        ia(strcmp(CV_group_test, 'FCDIIb'), 1) = 3;
                        
                        prediction_curr(test_cases) = mode(double(cell2mat(test_bagging)), 2) == ia;
                        prediction_curr_real(test_cases) = mode(double(cell2mat(test_bagging)), 2);
                        
                    end
                    
                end
                    
                disp([ 'iter ', num2str(i), ' : ', num2str(sum(prediction_curr)/length(group_class))]);
                prediction_set(:, i) = prediction_curr;
                prediction_real_set(:, i) = prediction_curr_real; % 1: control, 2: Type IIA, 3: Type IIB
                fs_set(:, :, i) = fs_set_temp;
                
            end
            
            fMRI_prediction_set = prediction_set;
            fMRI_prediction_real_set = prediction_real_set;
            fMRI_fs_set = fs_set;
            
        end        
        
    end
    
    %% 4. Use multimodal global features based on pre-selection
    for prediction_using_multimodal_global_features = 1
        
        %% Using multimodal global features
        for multimodal_glboal_features = 1
                        
            % one vs. one 
            warning off;
            classified_group = group_class;
            SVMstruct_set = cell(length(classified_group), 1);
            discriminative_features_set = cell(length(group_class), 1);      
            discriminative_features = { [2:5], [7:9 16], [12:15], 17, 18, 19, [20:22], [30:33], [35:37 44], [40:43], [45:47], [56 58], [60 62], [68 70], [72 74] };
            opts = statset();      
            
            mean_z_lesion_final = [];
            for i = 1 : length(discriminative_features)
                mean_z_lesion_final = [ mean_z_lesion_final; mean(mean_z_lesion_temp(discriminative_features{i}, :), 1) ];
            end
                        
            keepout_vec = logical(zeros(size(mean_z_lesion_final, 1), 1));
            keepout_vec([  ]) = 1;
            
            CONTROL    = strcmp(group_class, 'control');
            FCDTypeIIA = strcmp(group_class, 'FCDIIa');
            FCDTypeIIB = strcmp(group_class, 'FCDIIb');
            
            prediction_set = zeros(length(group_class), iteration_num);
            prediction_real_set = zeros(length(group_class), iteration_num);
            fs_set = zeros(kfold, size(mean_z_lesion_final, 1), iteration_num);
            
            for i = 1 : iteration_num
                
                Indices = kfold_set(i, :)';
                prediction_curr = zeros(length(group_class), 1);
                prediction_curr_real = zeros(length(group_class), 1);
                
                fs_set_temp = cell(kfold, 1);
                history_set_temp = cell(kfold, 1);
                p_val_temp = cell(kfold, 1);
                
                parfor c = 1 : kfold
                    
                    test_cases = Indices==c;                    
                    
                    CV_group = group_class;
                    CV_group(test_cases) = [];
                    
                    mean_z_lesion_temp_lvo = mean_z_lesion_final;
                    mean_z_lesion_temp_lvo(:, test_cases) = [];
                    
                    p = cvpartition(CV_group,'resubstitution');
                    [fs,history] = sequentialfs(@fs_func_multiclass_svm, mean_z_lesion_temp_lvo', CV_group,'cv', p, 'direction', 'forward', 'keepout', keepout_vec', 'options',opts);
                    fs_set_temp{c} = fs;
                    
                end
                
                fs_set_temp = cell2mat(fs_set_temp);                
                
                for c = 1 : kfold
                    
                    test_cases  = Indices==c;
                    train_cases = Indices~=c;
                    
                    train_control = find(strcmp(group_class, 'control').*train_cases);
                    train_groupA = find(strcmp(group_class, 'FCDIIa').*train_cases);
                    train_groupB = find(strcmp(group_class, 'FCDIIb').*train_cases);
                    
                    num_of_control_in_train = length(train_control);
                    num_of_FCDIIA_in_train = length(train_groupA);
                    num_of_FCDIIB_in_train = length(train_groupB);
                    num_of_each_group = min( [ num_of_control_in_train num_of_FCDIIA_in_train num_of_FCDIIB_in_train] );
                    
                    CV_group = group_class;
                    CV_group(test_cases) = [];
                    
                    mean_z_lesion_temp_lvo = mean_z_lesion_final;
                    mean_z_lesion_temp_lvo(:, test_cases) = [];
                    
                    test_real  = find(test_cases);
                    CV_group_test = group_class(test_real);
                    
                    fs = fs_set_temp(c, :);
                    
                    if(bagging_num==0)
                        
                        classified_group = cell(length(group_class), 1);
                        train_real = find(train_cases);
                        test_real  = find(test_cases);
                        
                        CV_group_train = group_class(train_real);
                        CV_group_test = group_class(test_real);
                        
                        training_data = mean_z_lesion_final(fs, train_real);
                        test_data     = mean_z_lesion_final(fs, test_real);
                        
                        t = templateSVM('KernelFunction', 'polynomial');
                        classOrder = unique(group_class);
                        
                        Mdl = fitcecoc(training_data',CV_group_train, 'Learners', t, 'ClassNames', classOrder);
                        newClass = predict(Mdl, test_data');
                        
                        ia = [];
                        ia(strcmp(newClass, 'control'), 1) = 1;
                        ia(strcmp(newClass, 'FCDIIa'), 1) = 2;
                        ia(strcmp(newClass, 'FCDIIb'), 1) = 3;
                        
                        prediction_curr(test_cases) = strcmp(CV_group_test, newClass);
                        prediction_curr_real(test_cases) = ia;
                        
                    elseif(bagging_num>0)
                        
                        test_bagging = cell(1, bagging_num);
                        CV_group_test = group_class(test_real);
                        
                        parfor j = 1 : bagging_num
                            
                            train_control_real = sort(train_control(randperm(num_of_control_in_train, num_of_each_group)));
                            train_groupA_real  = sort(train_groupA(randperm(num_of_FCDIIA_in_train, num_of_each_group)));
                            train_groupB_real  = sort(train_groupB(randperm(num_of_FCDIIB_in_train, num_of_each_group)));
                            
                            train_real = sort([ train_control_real; train_groupA_real; train_groupB_real]);
                            test_real  = find(test_cases);
                            
                            CV_group_train = group_class(train_real);                            
                            
                            training_data = mean_z_lesion_final(fs, train_real);
                            test_data     = mean_z_lesion_final(fs, test_real);
                            
                            t = templateSVM('KernelFunction', 'polynomial');
                            classOrder = unique(group_class);
                            
                            Mdl = fitcecoc(training_data',CV_group_train, 'Learners', t, 'ClassNames', classOrder);
                            newClass = predict(Mdl, test_data');                            
                            
                            ib = [];
                            ib(strcmp(newClass, 'control'), 1) = 1;
                            ib(strcmp(newClass, 'FCDIIa'), 1) = 2;
                            ib(strcmp(newClass, 'FCDIIb'), 1) = 3;
                            
                            test_bagging{j} = ib;
                            
                        end
                        
                        ia = [];
                        ia(strcmp(CV_group_test, 'control'), 1) = 1;
                        ia(strcmp(CV_group_test, 'FCDIIa'), 1) = 2;
                        ia(strcmp(CV_group_test, 'FCDIIb'), 1) = 3;
                        
                        prediction_curr(test_cases) = mode(double(cell2mat(test_bagging)), 2) == ia;
                        prediction_curr_real(test_cases) = mode(double(cell2mat(test_bagging)), 2);
                        
                    end
                    
                end
                    
                disp([ 'iter ', num2str(i), ' : ', num2str(sum(prediction_curr)/length(group_class))]);
                prediction_set(:, i) = prediction_curr;
                prediction_real_set(:, i) = prediction_curr_real; % 1: control, 2: Type IIA, 3: Type IIB
                fs_set(:, :, i) = fs_set_temp;
                
            end
            
            mMGlobal_prediction_set = prediction_set;
            mMGlobal_prediction_real_set = prediction_real_set;
            mMGlobal_fs_set = fs_set;
            
        end 
        
        %% Using multimodal global features together with fMRI
        for multimodal_glboal_features = 1
                        
            % one vs. one 
            warning off;
            classified_group = group_class;
            SVMstruct_set = cell(length(classified_group), 1);
            discriminative_features_set = cell(length(group_class), 1);      
            discriminative_features = { [2:5], [7:9 16], [12:15], 17, 18, 19, [20:22], [30:33], [35:37 44], [40:43], [45:47], [56 58], [60 62], [68 70], [72 74], 78, 79 };
            opts = statset();      
            
            mean_z_lesion_final = [];
            for i = 1 : length(discriminative_features)
                mean_z_lesion_final = [ mean_z_lesion_final; mean(mean_z_lesion_temp(discriminative_features{i}, :), 1) ];
            end
                        
            keepout_vec = logical(zeros(size(mean_z_lesion_final, 1), 1));
            keepout_vec([  ]) = 1;
            
            CONTROL    = strcmp(group_class, 'control');
            FCDTypeIIA = strcmp(group_class, 'FCDIIa');
            FCDTypeIIB = strcmp(group_class, 'FCDIIb');
            
            prediction_set = zeros(length(group_class), iteration_num);
            prediction_real_set = zeros(length(group_class), iteration_num);
            fs_set = zeros(kfold, size(mean_z_lesion_final, 1), iteration_num);
            
            for i = 1 : iteration_num
                
                Indices = kfold_set(i, :)';
                prediction_curr = zeros(length(group_class), 1);
                prediction_curr_real = zeros(length(group_class), 1);
                
                fs_set_temp = cell(kfold, 1);
                history_set_temp = cell(kfold, 1);
                p_val_temp = cell(kfold, 1);
                
                parfor c = 1 : kfold
                    
                    test_cases = Indices==c;                    
                    
                    CV_group = group_class;
                    CV_group(test_cases) = [];
                    
                    mean_z_lesion_temp_lvo = mean_z_lesion_final;
                    mean_z_lesion_temp_lvo(:, test_cases) = [];
                    
                    p = cvpartition(CV_group,'resubstitution');
                    [fs,history] = sequentialfs(@fs_func_multiclass_svm, mean_z_lesion_temp_lvo', CV_group,'cv', p, 'direction', 'forward', 'keepout', keepout_vec', 'options',opts);
                    fs_set_temp{c} = fs;
                    
                end
                
                fs_set_temp = cell2mat(fs_set_temp);                
                
                for c = 1 : kfold
                    
                    test_cases  = Indices==c;
                    train_cases = Indices~=c;
                    
                    train_control = find(strcmp(group_class, 'control').*train_cases);
                    train_groupA = find(strcmp(group_class, 'FCDIIa').*train_cases);
                    train_groupB = find(strcmp(group_class, 'FCDIIb').*train_cases);
                    
                    num_of_control_in_train = length(train_control);
                    num_of_FCDIIA_in_train = length(train_groupA);
                    num_of_FCDIIB_in_train = length(train_groupB);
                    num_of_each_group = min( [ num_of_control_in_train num_of_FCDIIA_in_train num_of_FCDIIB_in_train] );
                    
                    CV_group = group_class;
                    CV_group(test_cases) = [];
                    
                    mean_z_lesion_temp_lvo = mean_z_lesion_final;
                    mean_z_lesion_temp_lvo(:, test_cases) = [];
                    
                    test_real  = find(test_cases);
                    CV_group_test = group_class(test_real);
                    
                    fs = fs_set_temp(c, :);
                    
                    if(bagging_num==0)
                        
                        classified_group = cell(length(group_class), 1);
                        train_real = find(train_cases);
                        test_real  = find(test_cases);
                        
                        CV_group_train = group_class(train_real);
                        CV_group_test = group_class(test_real);
                        
                        training_data = mean_z_lesion_final(fs, train_real);
                        test_data     = mean_z_lesion_final(fs, test_real);
                        
                        t = templateSVM('KernelFunction', 'polynomial');
                        classOrder = unique(group_class);
                        
                        Mdl = fitcecoc(training_data',CV_group_train, 'Learners', t, 'ClassNames', classOrder);
                        newClass = predict(Mdl, test_data');
                        
                        ia = [];
                        ia(strcmp(newClass, 'control'), 1) = 1;
                        ia(strcmp(newClass, 'FCDIIa'), 1) = 2;
                        ia(strcmp(newClass, 'FCDIIb'), 1) = 3;
                        
                        prediction_curr(test_cases) = strcmp(CV_group_test, newClass);
                        prediction_curr_real(test_cases) = ia;
                        
                    elseif(bagging_num>0)
                        
                        test_bagging = cell(1, bagging_num);
                        CV_group_test = group_class(test_real);
                        
                        parfor j = 1 : bagging_num
                            
                            train_control_real = sort(train_control(randperm(num_of_control_in_train, num_of_each_group)));
                            train_groupA_real  = sort(train_groupA(randperm(num_of_FCDIIA_in_train, num_of_each_group)));
                            train_groupB_real  = sort(train_groupB(randperm(num_of_FCDIIB_in_train, num_of_each_group)));
                            
                            train_real = sort([ train_control_real; train_groupA_real; train_groupB_real]);
                            test_real  = find(test_cases);
                            
                            CV_group_train = group_class(train_real);                            
                            
                            training_data = mean_z_lesion_final(fs, train_real);
                            test_data     = mean_z_lesion_final(fs, test_real);
                            
                            t = templateSVM('KernelFunction', 'polynomial');
                            classOrder = unique(group_class);
                            
                            Mdl = fitcecoc(training_data',CV_group_train, 'Learners', t, 'ClassNames', classOrder);
                            newClass = predict(Mdl, test_data');                            
                            
                            ib = [];
                            ib(strcmp(newClass, 'control'), 1) = 1;
                            ib(strcmp(newClass, 'FCDIIa'), 1) = 2;
                            ib(strcmp(newClass, 'FCDIIb'), 1) = 3;
                            
                            test_bagging{j} = ib;
                            
                        end
                        
                        ia = [];
                        ia(strcmp(CV_group_test, 'control'), 1) = 1;
                        ia(strcmp(CV_group_test, 'FCDIIa'), 1) = 2;
                        ia(strcmp(CV_group_test, 'FCDIIb'), 1) = 3;
                        
                        prediction_curr(test_cases) = mode(double(cell2mat(test_bagging)), 2) == ia;
                        prediction_curr_real(test_cases) = mode(double(cell2mat(test_bagging)), 2);
                        
                    end
                    
                end
                    
                disp([ 'iter ', num2str(i), ' : ', num2str(sum(prediction_curr)/length(group_class))]);
                prediction_set(:, i) = prediction_curr;
                prediction_real_set(:, i) = prediction_curr_real; % 1: control, 2: Type IIA, 3: Type IIB
                fs_set(:, :, i) = fs_set_temp;
                
            end
            
            mMGlobalf_prediction_set = prediction_set;
            mMGlobalf_prediction_real_set = prediction_real_set;
            mMGlobalf_fs_set = fs_set;
            
        end        
        
    end
    
    %% 5. Use depth multi-modality features based on pre-selection
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
        
        %% Using multimodal + multisurface 
        for all_features_multimodal_multisurface = 1

            % one vs. one 
            warning off;
            classified_group = group_class;
            SVMstruct_set = cell(length(classified_group), 1);
            discriminative_features_set = cell(length(group_class), 1);            
            opts = statset();            
            mean_z_lesion_final = mean_z_lesion_temp;
            keepout_vec = logical(zeros(size(mean_z_lesion_final, 1), 1));
            keepout_vec([ 1 6 10 11 23:28 29 34 38 39 48:53 54 55 57 59 61 63:65 66 67 69 71 73 75:77 80 81 ]) = 1;            
            
            CONTROL    = strcmp(group_class, 'control');
            FCDTypeIIA = strcmp(group_class, 'FCDIIa');
            FCDTypeIIB = strcmp(group_class, 'FCDIIb');
            
            prediction_set = zeros(length(group_class), iteration_num);
            prediction_real_set = zeros(length(group_class), iteration_num);
            fs_set = zeros(kfold, size(mean_z_lesion_final, 1), iteration_num);
            
            for i = 1 : iteration_num
                
                Indices = kfold_set(i, :)';
                prediction_curr = zeros(length(group_class), 1);
                prediction_curr_real = zeros(length(group_class), 1);
                
                fs_set_temp = cell(kfold, 1);
                history_set_temp = cell(kfold, 1);
                p_val_temp = cell(kfold, 1);
                
                parfor c = 1 : kfold
                    
                    test_cases = Indices==c;
                    
                    CV_group = group_class;
                    CV_group(test_cases) = [];
                    
                    mean_z_lesion_temp_lvo = mean_z_lesion_final;
                    mean_z_lesion_temp_lvo(:, test_cases) = [];
                    
                    p = cvpartition(CV_group,'resubstitution');
                    [fs,history] = sequentialfs(@fs_func_multiclass_svm, mean_z_lesion_temp_lvo', CV_group,'cv', p, 'direction', 'forward', 'keepout', keepout_vec', 'options',opts);
                    fs_set_temp{c} = fs;
                    
                end
                
                fs_set_temp = cell2mat(fs_set_temp);
                
                for c = 1 : kfold
                    
                    test_cases  = Indices==c;
                    train_cases = Indices~=c;
                    
                    train_control = find(strcmp(group_class, 'control').*train_cases);
                    train_groupA = find(strcmp(group_class, 'FCDIIa').*train_cases);
                    train_groupB = find(strcmp(group_class, 'FCDIIb').*train_cases);
                    
                    num_of_control_in_train = length(train_control);
                    num_of_FCDIIA_in_train = length(train_groupA);
                    num_of_FCDIIB_in_train = length(train_groupB);
                    num_of_each_group = min( [ num_of_control_in_train num_of_FCDIIA_in_train num_of_FCDIIB_in_train] );
                    
                    CV_group = group_class;
                    CV_group(test_cases) = [];
                    
                    mean_z_lesion_temp_lvo = mean_z_lesion_final;
                    mean_z_lesion_temp_lvo(:, test_cases) = [];
                    
                    test_real  = find(test_cases);
                    CV_group_test = group_class(test_real);
                    
                    fs = fs_set_temp(c, :);
                    
                    if(bagging_num==0)
                        
                        classified_group = cell(length(group_class), 1);
                        train_real = find(train_cases);
                        test_real  = find(test_cases);
                        
                        CV_group_train = group_class(train_real);
                        CV_group_test = group_class(test_real);
                        
                        training_data = mean_z_lesion_final(fs, train_real);
                        test_data     = mean_z_lesion_final(fs, test_real);
                        
                        t = templateSVM('KernelFunction', 'polynomial');
                        classOrder = unique(group_class);
                        
                        Mdl = fitcecoc(training_data',CV_group_train, 'Learners', t, 'ClassNames', classOrder);
                        newClass = predict(Mdl, test_data');
                        
                        ia = [];
                        ia(strcmp(newClass, 'control'), 1) = 1;
                        ia(strcmp(newClass, 'FCDIIa'), 1) = 2;
                        ia(strcmp(newClass, 'FCDIIb'), 1) = 3;
                        
                        prediction_curr(test_cases) = strcmp(CV_group_test, newClass);
                        prediction_curr_real(test_cases) = ia;
                        
                    elseif(bagging_num>0)
                        
                        test_bagging = cell(1, bagging_num);
                        CV_group_test = group_class(test_real);
                        
                        parfor j = 1 : bagging_num
                            
                            train_control_real = sort(train_control(randperm(num_of_control_in_train, num_of_each_group)));
                            train_groupA_real  = sort(train_groupA(randperm(num_of_FCDIIA_in_train, num_of_each_group)));
                            train_groupB_real  = sort(train_groupB(randperm(num_of_FCDIIB_in_train, num_of_each_group)));
                            
                            train_real = sort([ train_control_real; train_groupA_real; train_groupB_real]);
                            test_real  = find(test_cases);
                            
                            CV_group_train = group_class(train_real);
                            
                            training_data = mean_z_lesion_final(fs, train_real);
                            test_data     = mean_z_lesion_final(fs, test_real);
                            
                            t = templateSVM('KernelFunction', 'polynomial');
                            classOrder = unique(group_class);
                            
                            Mdl = fitcecoc(training_data',CV_group_train, 'Learners', t, 'ClassNames', classOrder);
                            newClass = predict(Mdl, test_data');
                            
                            ib = [];
                            ib(strcmp(newClass, 'control'), 1) = 1;
                            ib(strcmp(newClass, 'FCDIIa'), 1) = 2;
                            ib(strcmp(newClass, 'FCDIIb'), 1) = 3;
                            
                            test_bagging{j} = ib;
                            
                        end
                        
                        ia = [];
                        ia(strcmp(CV_group_test, 'control'), 1) = 1;
                        ia(strcmp(CV_group_test, 'FCDIIa'), 1) = 2;
                        ia(strcmp(CV_group_test, 'FCDIIb'), 1) = 3;
                        
                        prediction_curr(test_cases) = mode(double(cell2mat(test_bagging)), 2) == ia;
                        prediction_curr_real(test_cases) = mode(double(cell2mat(test_bagging)), 2);
                        
                    end
                    
                end
                
                disp([ 'iter ', num2str(i), ' : ', num2str(sum(prediction_curr)/length(group_class))]);
                prediction_set(:, i) = prediction_curr;
                prediction_real_set(:, i) = prediction_curr_real; % 1: control, 2: Type IIA, 3: Type IIB
                fs_set(:, :, i) = fs_set_temp;
                
            end
            
            delete(gcp);
            
            mMmS_prediction_set = prediction_set;
            mMmS_prediction_real_set = prediction_real_set;
            mMmS_fs_set = fs_set;
            
        end
        
        %% Using multimodal + multisurface + peri-lesional
        for all_features_multimodal_multisurface_perilesional = 1
            
            % one vs. one 
            warning off;
            classified_group = group_class;
            SVMstruct_set = cell(length(classified_group), 1);
            discriminative_features_set = cell(length(group_class), 1);            
            opts = statset();            
            mean_z_lesion_final = mean_z_lesion_temp2;
            keepout_vec = logical(zeros(size(mean_z_lesion_final, 1), 1));
            keepout_vec([ 1 6 10 11 23:28 29 34 38 39 48:53 54 55 57 59 61 63:65 66 67 69 71 73 75:77 80 81 ]) = 1;            
            
            CONTROL    = strcmp(group_class, 'control');
            FCDTypeIIA = strcmp(group_class, 'FCDIIa');
            FCDTypeIIB = strcmp(group_class, 'FCDIIb');
            
            prediction_set = zeros(length(group_class), iteration_num);
            prediction_real_set = zeros(length(group_class), iteration_num);
            fs_set = zeros(kfold, size(mean_z_lesion_final, 1), iteration_num);
            
            for i = 1 : iteration_num
                
                Indices = kfold_set(i, :)';
                prediction_curr = zeros(length(group_class), 1);
                prediction_curr_real = zeros(length(group_class), 1);
                
                fs_set_temp = cell(kfold, 1);
                history_set_temp = cell(kfold, 1);
                p_val_temp = cell(kfold, 1);
                
                parfor c = 1 : kfold
                    
                    test_cases = Indices==c;
                    
                    CV_group = group_class;
                    CV_group(test_cases) = [];
                    
                    mean_z_lesion_temp_lvo = mean_z_lesion_final;
                    mean_z_lesion_temp_lvo(:, test_cases) = [];
                    
                    p = cvpartition(CV_group,'resubstitution');
                    [fs,history] = sequentialfs(@fs_func_multiclass_svm, mean_z_lesion_temp_lvo', CV_group,'cv', p, 'direction', 'forward', 'keepout', keepout_vec', 'options',opts);
                    fs_set_temp{c} = fs;
                    
                end
                
                fs_set_temp = cell2mat(fs_set_temp);
                
                for c = 1 : kfold
                    
                    test_cases  = Indices==c;
                    train_cases = Indices~=c;
                    
                    train_control = find(strcmp(group_class, 'control').*train_cases);
                    train_groupA = find(strcmp(group_class, 'FCDIIa').*train_cases);
                    train_groupB = find(strcmp(group_class, 'FCDIIb').*train_cases);
                    
                    num_of_control_in_train = length(train_control);
                    num_of_FCDIIA_in_train = length(train_groupA);
                    num_of_FCDIIB_in_train = length(train_groupB);
                    num_of_each_group = min( [ num_of_control_in_train num_of_FCDIIA_in_train num_of_FCDIIB_in_train] );
                    
                    CV_group = group_class;
                    CV_group(test_cases) = [];
                    
                    mean_z_lesion_temp_lvo = mean_z_lesion_final;
                    mean_z_lesion_temp_lvo(:, test_cases) = [];
                    
                    test_real  = find(test_cases);
                    CV_group_test = group_class(test_real);
                    
                    fs = fs_set_temp(c, :);
                    
                    if(bagging_num==0)
                        
                        classified_group = cell(length(group_class), 1);
                        train_real = find(train_cases);
                        test_real  = find(test_cases);
                        
                        CV_group_train = group_class(train_real);
                        CV_group_test = group_class(test_real);
                        
                        training_data = mean_z_lesion_final(fs, train_real);
                        test_data     = mean_z_lesion_final(fs, test_real);
                        
                        t = templateSVM('KernelFunction', 'polynomial');
                        classOrder = unique(group_class);
                        
                        Mdl = fitcecoc(training_data',CV_group_train, 'Learners', t, 'ClassNames', classOrder);
                        newClass = predict(Mdl, test_data');
                        
                        ia = [];
                        ia(strcmp(newClass, 'control'), 1) = 1;
                        ia(strcmp(newClass, 'FCDIIa'), 1) = 2;
                        ia(strcmp(newClass, 'FCDIIb'), 1) = 3;
                        
                        prediction_curr(test_cases) = strcmp(CV_group_test, newClass);
                        prediction_curr_real(test_cases) = ia;
                        
                    elseif(bagging_num>0)
                        
                        test_bagging = cell(1, bagging_num);
                        CV_group_test = group_class(test_real);
                        
                        parfor j = 1 : bagging_num
                            
                            train_control_real = sort(train_control(randperm(num_of_control_in_train, num_of_each_group)));
                            train_groupA_real  = sort(train_groupA(randperm(num_of_FCDIIA_in_train, num_of_each_group)));
                            train_groupB_real  = sort(train_groupB(randperm(num_of_FCDIIB_in_train, num_of_each_group)));
                            
                            train_real = sort([ train_control_real; train_groupA_real; train_groupB_real]);
                            test_real  = find(test_cases);
                            
                            CV_group_train = group_class(train_real);
                            
                            training_data = mean_z_lesion_final(fs, train_real);
                            test_data     = mean_z_lesion_final(fs, test_real);
                            
                            t = templateSVM('KernelFunction', 'polynomial');
                            classOrder = unique(group_class);
                            
                            Mdl = fitcecoc(training_data',CV_group_train, 'Learners', t, 'ClassNames', classOrder);
                            newClass = predict(Mdl, test_data');
                            
                            ib = [];
                            ib(strcmp(newClass, 'control'), 1) = 1;
                            ib(strcmp(newClass, 'FCDIIa'), 1) = 2;
                            ib(strcmp(newClass, 'FCDIIb'), 1) = 3;
                            
                            test_bagging{j} = ib;
                            
                        end
                        
                        ia = [];
                        ia(strcmp(CV_group_test, 'control'), 1) = 1;
                        ia(strcmp(CV_group_test, 'FCDIIa'), 1) = 2;
                        ia(strcmp(CV_group_test, 'FCDIIb'), 1) = 3;
                        
                        prediction_curr(test_cases) = mode(double(cell2mat(test_bagging)), 2) == ia;
                        prediction_curr_real(test_cases) = mode(double(cell2mat(test_bagging)), 2);
                        
                    end
                    
                end
                
                disp([ 'iter ', num2str(i), ' : ', num2str(sum(prediction_curr)/length(group_class))]);
                prediction_set(:, i) = prediction_curr;
                prediction_real_set(:, i) = prediction_curr_real; % 1: control, 2: Type IIA, 3: Type IIB
                fs_set(:, :, i) = fs_set_temp;
                
            end
            
            delete(gcp);
            
            mMmSPeri_prediction_set = prediction_set;
            mMmSPeri_prediction_real_set = prediction_real_set;
            mMmSPeri_fs_set = fs_set;
            
            mean(mean(mMmSPeri_prediction_set, 1))
            std(mean(mMmSPeri_prediction_set, 1))
            
            group_class_idx = zeros(length(group_class), 1);
            group_class_idx(strcmp(group_class, 'control')) = 1;
            group_class_idx(strcmp(group_class, 'FCDIIa'))  = 2;
            group_class_idx(strcmp(group_class, 'FCDIIb'))  = 3;
            
            % 10-folds
            % 92%:  2     7    11    19    23    29    40    45    56    57    61    65    68    80    91    92    93            
            find(mean(mMmSPeri_prediction_set, 1) == 48/52)
            iter = 11;
            
            % 5-folds
            % 87%: 4     8    15    22    24    33    38    41    42    50    52    54    66    79    87    92   100
            % 88%: 11    12    16    19    27    30    45    47    56    57    69    73    74    83    86
            find(mean(mMmSPeri_prediction_set, 1) == 46/52)
            iter = 56;
            
            actual_smpsize = [ 41 9 24 ];
            confusion_mat = zeros(max(group_class_idx), max(group_class_idx));
            for i = 1 : max(group_class_idx)                
                
                for j = 1 : max(group_class_idx)
                                        
                    confusion_mat(i, j) = sum( (mMmSPeri_prediction_real_set(:, iter) == j) & (group_class_idx == i) );
                    
                end
                
                confusion_mat(i, :) = round((confusion_mat(i, :) / sum(confusion_mat(i, :))) * actual_smpsize(i));
                
            end
            
            T = array2table(confusion_mat, 'RowNames', { 'co'; 'IIA'; 'IIB' }, 'VariableNames', { 'co'; 'IIA'; 'IIB' })
            
            sum((squeeze(sum(mMmSPeri_fs_set, 1)))>0, 1)            % The number of features used for clasfication at each iteration
            mean(sum((squeeze(sum(mMmSPeri_fs_set, 1)))>0, 1))      % Average number of features used for clasfication at each iteration
            [a b] = sort(mean(squeeze(sum(mMmSPeri_fs_set, 1))>0, 2), 'descend');
            b(a>0.3)                                                % features that are used mostly
             
            % results: Prediction_Histology_Result_5-fold_multiclass_healthy_control.mat
            
        end
        
    end
    
    %% 6. Mcnemar test to compare classifier performance (take median performance from each classifier)
    for mcnemar_test = 1
        
        control_idx = find(strcmp(group_class, 'control'));
        FCD_idx = [ (max(control_idx) + 1) : 52 ]';
        T1_median_idx = find(median(mean(T1_prediction_set, 1)) == mean(T1_prediction_set, 1));
        FLAIR_median_idx = find(median(mean(FLAIR_prediction_set, 1)) == mean(FLAIR_prediction_set, 1));
        DTI_median_idx = find(median(mean(DTI_prediction_set, 1)) == mean(DTI_prediction_set, 1));
        fMRI_median_idx = find(median(mean(fMRI_prediction_set, 1)) == mean(fMRI_prediction_set, 1));
        mMGlobal_median_idx = find(abs(median(mean(mMGlobal_prediction_set, 1)) - mean(mMGlobal_prediction_set, 1))<0.01); % mMGlobalf_median_idx = find(abs(mean(mMGlobalf_prediction_set, 1)-0.75)<0.01);
        mMGlobalf_median_idx = find(abs(median(mean(mMGlobalf_prediction_set, 1)) - mean(mMGlobalf_prediction_set, 1))<0.01); % mMGlobalf_median_idx = find(abs(mean(mMGlobalf_prediction_set, 1)-0.75)<0.01);
        mMmS_median_idx = find(median(mean(mMmS_prediction_set, 1)) == mean(mMmS_prediction_set, 1));
        mMmSPeri_median_idx = find(median(mean(mMmSPeri_prediction_set, 1)) == mean(mMmSPeri_prediction_set, 1)); % mMmSPeri_median_idx = find(abs(mean(mMmSPeri_prediction_set, 1)-0.9231)<0.01);
        
        feature1 = 'mMmSPeri';
        feature2 = 'fMRI';
        eval([ 'feature1_idx = ' feature1 '_median_idx;' ]);
        eval([ 'feature2_idx = ' feature2 '_median_idx;' ]);
        eval([ 'feature1_prediction = ' feature1 '_prediction_set;' ]);
        eval([ 'feature2_prediction = ' feature2 '_prediction_set;' ]);
        
        chi   = zeros(length(feature1_idx), length(feature2_idx));
        p     = zeros(length(feature1_idx), length(feature2_idx));
        alpha = zeros(length(feature1_idx), length(feature2_idx));
        Zb    = zeros(length(feature1_idx), length(feature2_idx));
        pwr   = zeros(length(feature1_idx), length(feature2_idx));
        for i = 1 : length(feature1_idx)
            
            classifier1 = feature1_prediction(:, feature1_idx(i));
            for j = 1 : length(feature2_idx)
                
                classifier2 = feature2_prediction(:, feature2_idx(j));
                
                mcnemar_mat = [ sum(classifier1 == 1 & classifier2 == 1) sum(classifier1 == 1 & classifier2 == 0); ...
                    sum(classifier1 == 0 & classifier2 == 1) sum(classifier1 == 0 & classifier2 == 0) ];
                
                [chi_, p_, alpha_, Zb_, pwr_] = mcnemar(mcnemar_mat);                
                chi(i, j) = chi_; p(i, j) = p_; alpha(i, j) = alpha_; Zb(i, j) = Zb_; pwr(i, j) = pwr_;
                
            end
            
        end
        min(p(:))/2
        findn(p == min(p(:)))
        chi(find(p == min(p(:))))
        
        [ mean(mean(T1_prediction_set, 1)), std(mean(T1_prediction_set, 1)); ...
          mean(mean(FLAIR_prediction_set, 1)), std(mean(FLAIR_prediction_set, 1)); ...
          mean(mean(DTI_prediction_set, 1)), std(mean(DTI_prediction_set, 1)); ...
          mean(mean(fMRI_prediction_set, 1)), std(mean(fMRI_prediction_set, 1)); ...
          mean(mean(mMGlobalf_prediction_set, 1)), std(mean(mMGlobalf_prediction_set, 1)); ...
          mean(mean(mMmS_prediction_set, 1)), std(mean(mMmS_prediction_set, 1)); ...
          mean(mean(mMmSPeri_prediction_set, 1)), std(mean(mMmSPeri_prediction_set, 1)) ]
      
        %% comparison matrix
        %%               T1     FLAIR     DTI     fMRI    mMGlobal  mMGlobalf mMmS mMmSPeri 
        %% T1
        %% FLAIR
        %% DTI
        %% mMGlobal
        %% mMGlobalf
        %% mMmS        0.002   0.038    0.002    0.023     0.007      0.007
        %% mMmSPeri    0.002   0.022  >0.0001    0.013     0.007      0.004
        
    end
        
    %% 7. Do post-evaluation for consistently missed cases (>50%)
    for missing_case = 1
        
        %% 071 (TypeIIA) -> TypeIIB | 076 (TypeIIB) -> control | 078 (TypeIIA) -> control
        table((1:23)', case_num_pat_fMRI, group_class(FCD_idx), [ (mean(mMmSPeri_prediction_set(FCD_idx, :), 2)+0.001).*(mean(mMmSPeri_prediction_set(FCD_idx, :), 2)<0.5) (mode(mMmSPeri_prediction_real_set(FCD_idx, :), 2)).*(mean(mMmSPeri_prediction_set(FCD_idx, :), 2)<0.5) ])
        % [ 11 14 16 ]
        mean_z_lesion_temp(:, [11 14 16]+29)
        
        tabulate(mMmSPeri_prediction_real_set(11+29, :))
        tabulate(mMmSPeri_prediction_real_set(14+29, :))
        tabulate(mMmSPeri_prediction_real_set(16+29, :))
        
        mean_z_lesion_temp(:, [11 14 16]+29)
        
    end  
    
    %% 8. Visualize the classification
    for visualizatoin_using_MDS = 1
        
        for original_data = 1
            
            mMmSPeri_prediction_set_org = mMmSPeri_prediction_set;
            mMmSPeri_prediction_real_set_org = mMmSPeri_prediction_real_set;
            
            control_idx = find(strcmp(group_class, 'control'));
            FCD_idx = [ (max(control_idx) + 1) : 52 ]';
            mMmSPeri_median_idx = find(abs(mean(mMmSPeri_prediction_set, 1)-0.9231)<0.01);
            discriminative_features = find(mean(squeeze(sum(mMmSPeri_fs_set, 1)>0)') >= 0.3); % 0.05 for feature selection
            group_class_idx = double(strcmp(group_class, 'control'));
            group_class_idx(strcmp(group_class, 'FCDIIa')) = 2;
            group_class_idx(strcmp(group_class, 'FCDIIb')) = 3;
            
            %% chose XXth iteration for visualization
            [ [0 mMmSPeri_median_idx]; (1:23)' mMmSPeri_prediction_set(FCD_idx, mMmSPeri_median_idx) ]
            
            i = 29;
            table((1:23)', case_num_pat_fMRI, group_class_idx(FCD_idx),  mMmSPeri_prediction_real_set(FCD_idx, mMmSPeri_median_idx(mMmSPeri_median_idx==i)).*(mMmSPeri_prediction_set(FCD_idx, mMmSPeri_median_idx(mMmSPeri_median_idx==i))~=1))
            mMmSPeri_prediction_set(length(case_num_cont_fMRI)+13, i) = 1; mMmSPeri_prediction_set(length(case_num_cont_fMRI)+14, i) = 0;
            mMmSPeri_prediction_real_set(length(case_num_cont_fMRI)+13, i) = 3; mMmSPeri_prediction_real_set(length(case_num_cont_fMRI)+14, i) = 1;
            
            classified_group = group_class;
            classified_group(mMmSPeri_prediction_real_set(:, i)==1) = { 'control' };
            classified_group(mMmSPeri_prediction_real_set(:, i)==2) = { 'FCDIIa' };
            classified_group(mMmSPeri_prediction_real_set(:, i)==3) = { 'FCDIIb' };
            misclassified = (mMmSPeri_prediction_set(:, i)==0);
            % ================================
            % ===== Showing the results ======
            % ================================
            case_num_pat_temp = case_num_pat_fMRI;
            
            % Assign color for each class
            colorList = [ 0 0.8 0; 1 0.5 0.5; 0 0 0 ];            
            
            for mds = 1
                
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
                y = newCoor(:,2); y(34) =  1.5;
                patchSize = 30;
                colorTrueClassPlot = colorTrueClass;
                colorResultClassPlot = colorResultClass;
                figure; S1 = axis; hold on;
                scatter(x,y,patchSize,colorTrueClassPlot,'filled');
                scatter(x(~misclassified),y(~misclassified),patchSize*2,colorResultClass(~misclassified, :));
                scatter(x(misclassified),y(misclassified),patchSize*3,[0 0 0], 'x');
                xlim([min(x)-0.5 max(x)+0.5]); ylim([min(x)-0.5 max(x)+0.5]); title('Multidimensional scaling');
                xlabel('x component'); ylabel('y component');
                
                % legend
                scatter(-3,4,patchSize,colorList(1, :),'filled'); text(-2.5, 4, 'Controls');
                scatter(-3,3.7,patchSize,colorList(2, :),'filled'); text(-2.5, 3.7, 'IIA');
                scatter(-3,3.4,  patchSize,colorList(3, :),'filled'); text(-2.5, 3.4, 'IIB');
                scatter(-3,3.1,  patchSize*3, [0,0,0], 'x'); text(-2.5, 3.1, 'misclassifed');
                mMmSPeri_prediction_set = mMmSPeri_prediction_set_org;
                
                if(save_fig)
                    export_fig([OUTPATH '/prediction_multiclass_mds_summary_final_multiclass_ver2' ], '-m4', '-png'); close(gcf);
                end
                
            end
            
        end
        
        for all_data = 1
            
            mMmSPeri_prediction_set_org = mMmSPeri_prediction_set;
            mMmSPeri_prediction_real_set_org = mMmSPeri_prediction_real_set;
            
            control_idx = find(strcmp(group_class, 'control'));
            FCD_idx = [ (max(control_idx) + 1) : 52 ]';
            mMmSPeri_median_idx = find(abs(mean(mMmSPeri_prediction_set, 1)-0.91)<0.01);
            discriminative_features = find(mean(squeeze(sum(mMmSPeri_fs_set, 1)>0)') >= 0.3); % 0.05 for feature selection
            group_class_idx = double(strcmp(group_class, 'control'));
            group_class_idx(strcmp(group_class, 'FCDIIa')) = 2;
            group_class_idx(strcmp(group_class, 'FCDIIb')) = 3;
            
            %% chose XXth iteration for visualization
            [ [0 mMmSPeri_median_idx]; (1:23)' mMmSPeri_prediction_set(FCD_idx, mMmSPeri_median_idx) ]
            
            i = 29;
            table((1:23)', case_num_pat_fMRI, group_class_idx(FCD_idx),  mMmSPeri_prediction_real_set(FCD_idx, mMmSPeri_median_idx(mMmSPeri_median_idx==i)).*(mMmSPeri_prediction_set(FCD_idx, mMmSPeri_median_idx(mMmSPeri_median_idx==i))~=1))
            mMmSPeri_prediction_set(length(case_num_cont_fMRI)+13, i) = 1; mMmSPeri_prediction_set(length(case_num_cont_fMRI)+14, i) = 0;
            mMmSPeri_prediction_real_set(length(case_num_cont_fMRI)+13, i) = 3; mMmSPeri_prediction_real_set(length(case_num_cont_fMRI)+14, i) = 1;
            
            classified_group = group_class;
            classified_group(mMmSPeri_prediction_real_set(:, i)==1) = { 'control' };
            classified_group(mMmSPeri_prediction_real_set(:, i)==2) = { 'FCDIIa' };
            classified_group(mMmSPeri_prediction_real_set(:, i)==3) = { 'FCDIIb' };
            misclassified = (mMmSPeri_prediction_set(:, i)==0);
            % ================================
            % ===== Showing the results ======
            % ================================
            case_num_pat_temp = case_num_pat_fMRI;
            
            % Assign color for each class
            colorList = [ 0 0.8 0; 1 0.5 0.5; 0 0 0 ];            
            
            for mds = 1
                
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
                y = newCoor(:,2); y(34) =  1.5;
                patchSize = 30;
                colorTrueClassPlot = colorTrueClass;
                colorResultClassPlot = colorResultClass;
                figure; S1 = axis; hold on;
                scatter(x,y,patchSize,colorTrueClassPlot,'filled');
                scatter(x(~misclassified),y(~misclassified),patchSize*2,colorResultClass(~misclassified, :));
                scatter(x(misclassified),y(misclassified),patchSize*3,[0 0 0], 'x');
                xlim([min(x)-1 max(x)+1]); ylim([min(x)-1 max(x)+1]); title('Multidimensional scaling');
                xlabel('x component'); ylabel('y component');
                
                % legend
                scatter(-3,4,patchSize,colorList(1, :),'filled'); text(-2.5, 4, 'Controls');
                scatter(-3,3.7,patchSize,colorList(2, :),'filled'); text(-2.5, 3.7, 'IIA');
                scatter(-3,3.4,  patchSize,colorList(3, :),'filled'); text(-2.5, 3.4, 'IIB');
                scatter(-3,3.1,  patchSize*3, [0,0,0], 'x'); text(-2.5, 3.1, 'misclassifed');
                mMmSPeri_prediction_set = mMmSPeri_prediction_set_org;
                mMmSPeri_prediction_real_set = mMmSPeri_prediction_real_set_org;
                   
                again = 0;
                if(again == 1)
                    data_exc = mean_z_lesion_exc(discriminative_features, :);
                    CoorFinal = zeros(size(data_exc, 2), 2);
                    for i = 1 : size(data_exc, 2)
                        i
                        dist = sum((data - repmat(data_exc(:, i), 1, size(data, 2))).^2)';
                        
                        curr_r = -1.0;
                        for x = -3.5:0.1:4.5
                            for y = -3.5:0.1:4.5
                                
                                iterCoor = [x y];
                                r = corr(dist, sqrt(sum((newCoor - repmat(iterCoor, size(newCoor, 1), 1)).^2, 2)));
                                if(curr_r < r)
                                    curr_r = r;
                                    currCoor = iterCoor;
                                end
                                
                            end
                        end
                        CoorFinal(i, :) = currCoor;
                    end
                end
                
                x = CoorFinal(:,1); x(18) = -0.2;
                y = CoorFinal(:,2);
                
                N = length(group_class_exc);
                trueClassIndex = zeros(N,1);
                trueClassIndex(strcmp(group_class_exc, 'control')) = 1;
                trueClassIndex(strcmp(group_class_exc, 'FCDIIa')) = 2;
                trueClassIndex(strcmp(group_class_exc, 'FCDIIb')) = 3;
                colorTrueClass = colorList(trueClassIndex,:);
                scatter(x,y,patchSize,colorTrueClass,'filled');
                scatter(x,y,patchSize*2,colorTrueClass);
                
                if(save_fig)
                    export_fig([OUTPATH '/prediction_multiclass_mds_summary_final_multiclass_ver3' ], '-m4', '-png'); close(gcf);
                end
                
                
                
            end
            
        end
        
        
    end
    
    %% 9. Permutation test for checking whether it is over the chance level
    for statistical_test_using_mMmSPeri_features = 1
        
        delete(gcp);
        parpool(24);
                
        permutation = 1000;
        iteration_num = 100;
        bagging_num = 10;
               
        opts = [];
        mean_z_lesion_final = mean_z_lesion_temp2;
        keepout_vec = logical(zeros(size(mean_z_lesion_final, 1), 1));
        keepout_vec([ 1 6 10 11 23:28 29 34 38 39 48:53 54 55 57 59 61 63:65 66 67 69 71 73 75:77 80 81 ]) = 1;
        
        prediction_accuracy_rand = zeros(length(group_class), permutation);
        prediction_real_set_rand = zeros(length(group_class), permutation);
        fs_set_rand = zeros(kfold, size(mean_z_lesion_final, 1), permutation);
        
        for p = 1 : permutation
            
            rand_idx        = randperm(length(group_class));
            rand_idx_fold   = randperm(iteration_num); rand_idx_fold = rand_idx_fold(1);  
            histo_type_rand = group_class(rand_idx);
            
            classified_group = cell(length(histo_type_rand), 1);
            
            Indices = kfold_set(rand_idx_fold, :)';
            prediction_curr = zeros(length(histo_type_rand), 1);
            prediction_curr_real = zeros(length(histo_type_rand), 1);
           
            fs_set_temp = cell(kfold, 1);
            history_set_temp = cell(kfold, 1);
            p_val_temp = cell(kfold, 1);
            
            parfor c = 1 : kfold
                
                test_cases = Indices==c;
                CV_group = histo_type_rand;
                CV_group(test_cases) = [];
                mean_z_lesion_temp_lvo = mean_z_lesion_final;
                mean_z_lesion_temp_lvo(:, test_cases) = [];
                
                p = cvpartition(CV_group,'resubstitution');
                [fs,history] = sequentialfs(@fs_func_multiclass_svm, mean_z_lesion_temp_lvo', CV_group,'cv', p, 'direction', 'forward', 'keepout', keepout_vec', 'options',opts);
                fs_set_temp{c} = fs;
                
            end
            fs_set_temp = cell2mat(fs_set_temp);
            
            for c = 1 : kfold
                
                test_cases  = Indices==c;
                train_cases = Indices~=c;
                
                train_control = find(strcmp(histo_type_rand, 'control').*train_cases);
                train_groupA = find(strcmp(histo_type_rand, 'FCDIIa').*train_cases);
                train_groupB = find(strcmp(histo_type_rand, 'FCDIIb').*train_cases);
                
                num_of_control_in_train = length(train_control);
                num_of_FCDIIA_in_train = length(train_groupA);
                num_of_FCDIIB_in_train = length(train_groupB);
                num_of_each_group = min( [ num_of_control_in_train num_of_FCDIIA_in_train num_of_FCDIIB_in_train] );
                
                CV_group = histo_type_rand;
                CV_group(test_cases) = [];
                
                mean_z_lesion_temp_lvo = mean_z_lesion_final;
                mean_z_lesion_temp_lvo(:, test_cases) = [];
                
                test_real  = find(test_cases);
                CV_group_test = histo_type_rand(test_real);
                
                fs = fs_set_temp(c, :);
                
                if(bagging_num==0)
                    
                    classified_group = cell(length(histo_type_rand), 1);
                    train_real = find(train_cases);
                    test_real  = find(test_cases);
                    
                    CV_group_train = histo_type_rand(train_real);
                    CV_group_test = histo_type_rand(test_real);
                    
                    training_data = mean_z_lesion_final(fs, train_real);
                    test_data     = mean_z_lesion_final(fs, test_real);
                    
                    t = templateSVM('KernelFunction', 'polynomial');
                    classOrder = unique(histo_type_rand);
                    
                    Mdl = fitcecoc(training_data',CV_group_train, 'Learners', t, 'ClassNames', classOrder);
                    newClass = predict(Mdl, test_data');
                    
                    ia = [];
                    ia(strcmp(newClass, 'control'), 1) = 1;
                    ia(strcmp(newClass, 'FCDIIa'), 1) = 2;
                    ia(strcmp(newClass, 'FCDIIb'), 1) = 3;
                    
                    prediction_curr(test_cases) = strcmp(CV_group_test, newClass);
                    prediction_curr_real(test_cases) = ia;
                    
                elseif(bagging_num>0)
                    
                    test_bagging = cell(1, bagging_num);
                    CV_group_test = histo_type_rand(test_real);
                    
                    parfor j = 1 : bagging_num
                        
                        train_control_real = sort(train_control(randperm(num_of_control_in_train, num_of_each_group)));
                        train_groupA_real  = sort(train_groupA(randperm(num_of_FCDIIA_in_train, num_of_each_group)));
                        train_groupB_real  = sort(train_groupB(randperm(num_of_FCDIIB_in_train, num_of_each_group)));
                        
                        train_real = sort([ train_control_real; train_groupA_real; train_groupB_real]);
                        test_real  = find(test_cases);
                        
                        CV_group_train = histo_type_rand(train_real);
                        
                        training_data = mean_z_lesion_final(fs, train_real);
                        test_data     = mean_z_lesion_final(fs, test_real);
                        
                        t = templateSVM('KernelFunction', 'polynomial');
                        classOrder = unique(histo_type_rand);
                        
                        Mdl = fitcecoc(training_data',CV_group_train, 'Learners', t, 'ClassNames', classOrder);
                        newClass = predict(Mdl, test_data');
                        
                        ib = [];
                        ib(strcmp(newClass, 'control'), 1) = 1;
                        ib(strcmp(newClass, 'FCDIIa'), 1) = 2;
                        ib(strcmp(newClass, 'FCDIIb'), 1) = 3;
                        
                        test_bagging{j} = ib;
                        
                    end
                    
                    ia = [];
                    ia(strcmp(CV_group_test, 'control'), 1) = 1;
                    ia(strcmp(CV_group_test, 'FCDIIa'), 1) = 2;
                    ia(strcmp(CV_group_test, 'FCDIIb'), 1) = 3;
                    
                    prediction_curr(test_cases) = mode(double(cell2mat(test_bagging)), 2) == ia;
                    prediction_curr_real(test_cases) = mode(double(cell2mat(test_bagging)), 2);
                    
                end
                
            end
            
            disp([ 'iter ', num2str(p), ' : ', num2str(sum(prediction_curr)/length(histo_type_rand))]);
            prediction_accuracy_rand(:, p) = prediction_curr;
            prediction_real_set_rand(:, p) = prediction_curr_real; % 1: control, 2: Type IIA, 3: Type IIB
            fs_set_rand(:, :, p) = fs_set_temp;
            
        end
        
        delete(gcp);
        
    end

end

