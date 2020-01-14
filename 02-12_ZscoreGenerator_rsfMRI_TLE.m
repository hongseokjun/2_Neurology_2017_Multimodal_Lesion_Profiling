clear all
close all
addpath('/data/noel/noel2/local/matlab/surfstat_chicago');

%% Read demograpic data
for read_demodata = 1
    
    OutDir                    = '/local_raid/seokjun/01_project/04_IntracorticalAnalysis/03_result/'
    
    Cases_pat                 = '/local_raid/seokjun/01_project/04_IntracorticalAnalysis/01_analysis/Demographic_data_FCD_fMRI.txt'
    Group_pat                 = 'FCD'
    Prefix_pat                = 'mcd'
    
    Cases_healthy_cont        = '/local_raid/seokjun/01_project/04_IntracorticalAnalysis/01_analysis/Demographic_data_control_fMRI.txt'    
    Group_healthy_cont        = 'control'
    Prefix_healthy_cont       = 'TLE'
    
    Cases_disease_cont        = '/local_raid/seokjun/01_project/04_IntracorticalAnalysis/01_analysis/Demographic_data_TLE_final.txt'
    Group_disease_cont        = 'TLE'
    Prefix_disease_cont       = 'TLE'
        
    Left_Right                = 'both'
    NumIntSurf                = 3
    NumSubSurf                = 7
    NumMesh                   = 81920
    SamplingSpace             = 'native'
    img_contrast              = {'t1'}
    Kernel                    = 5
    Parametric                = 'quadratic'
    average_surface_dir       = '/local_raid/seokjun/01_project/04_IntracorticalAnalysis/03_result/average_surfaces/'
    
    visualization = 1;
    load('/local_raid/seokjun/01_project/04_IntracorticalAnalysis/01_analysis/pipeline/postprocessing_zscore_visualization/colormap_noel1.mat');
    
    fid                         = fopen(Cases_pat);
    demo                        = textscan(fid, '%s%f%s%s%s', 'Delimiter', ',', 'CollectOutput', 1);
    case_num_pat                = demo{1};
    age_pat                     = demo{2};
    gender_pat                  = demo{3}(:, 1);
    histo_type                  = demo{3}(:, 3);
    fclose(fid);
    
    fid                         = fopen(Cases_healthy_cont);
    demo                        = textscan(fid, '%s%f%s%s%f', 'Delimiter', ',', 'CollectOutput', 1);
    case_num_healthy_cont       = demo{1};
    age_healthy_cont            = demo{2};
    gender_healthy_cont         = demo{3}(:, 1);
    fclose(fid);
    
    fid = fopen(Cases_disease_cont);
    demo = textscan(fid, '%s%f%s%s', 'Delimiter', ',', 'CollectOutput', 1);
    case_num_disease_cont       = demo{1};
    age_disease_cont            = demo{2};
    gender_disease_cont         = demo{3}(:, 1);
    seizure_lateralization  = demo{3}(:, 2);
    fclose(fid);
    
    mean_age_healthy_cont   = mean(age_healthy_cont);
    std_age_healthy_cont    = std(age_healthy_cont, 0);
    
    mean_age_disease_cont   = mean(age_disease_cont);
    std_age_disease_cont    = std(age_disease_cont, 0);
    
    [h,p,ci,stats] = ttest2(age_healthy_cont,age_disease_cont);
    
end

%% Set up some parameters and variables
for setup_param = 1
    
    if(NumMesh == 81920)
        vertexnum = 81924;
    elseif(NumMesh == 327680)
        vertexnum = 327684;
    end
    
    if(strcmp(Left_Right, 'left'))
        vertexnum = vertexnum/2;
        Left_Right = { 'left'; 1:vertexnum };
    elseif(strcmp(Left_Right, 'right'))
        vertexnum = vertexnum/2;
        Left_Right = { 'right'; vertexnum+1:(vertexnum*2) };
    elseif(strcmp(Left_Right, 'both'))
        Left_Right = {'left', 'right'; 1:vertexnum/2 vertexnum/2+1:vertexnum };
    end
    
    TemplateSurf = cell(NumSubSurf+NumIntSurf+2, 1);
    TemplateSurf{1} = SurfStatReadSurf({[ average_surface_dir 'average_gray_surface_left_',  num2str(NumMesh) '_t1.obj' ],  ...
        [ average_surface_dir 'average_gray_surface_right_', num2str(NumMesh) '_t1.obj' ] });
    TemplateSurf{2} = SurfStatReadSurf({[ average_surface_dir 'average_intracortical_surface_1_left_',  num2str(NumMesh) '_t1.obj' ],  ...
        [ average_surface_dir 'average_intracortical_surface_1_right_', num2str(NumMesh) '_t1.obj' ] });
    TemplateSurf{3} = SurfStatReadSurf({[ average_surface_dir 'average_intracortical_surface_2_left_',  num2str(NumMesh) '_t1.obj' ],  ...
        [ average_surface_dir 'average_intracortical_surface_2_right_', num2str(NumMesh) '_t1.obj' ] });
    TemplateSurf{4} = SurfStatReadSurf({[ average_surface_dir 'average_intracortical_surface_3_left_',  num2str(NumMesh) '_t1.obj' ],  ...
        [ average_surface_dir 'average_intracortical_surface_3_right_', num2str(NumMesh) '_t1.obj' ] });
    TemplateSurf{5} = SurfStatReadSurf({[ average_surface_dir 'average_white_surface_left_',  num2str(NumMesh) '_t1.obj' ],  ...
        [ average_surface_dir 'average_white_surface_right_', num2str(NumMesh) '_t1.obj' ] });
    TemplateSurf{6} = SurfStatReadSurf({[ average_surface_dir 'average_white_surface_1_left_',  num2str(NumMesh) '_t1.obj' ],  ...
        [ average_surface_dir 'average_white_surface_1_right_', num2str(NumMesh) '_t1.obj' ] });
    TemplateSurf{7} = SurfStatReadSurf({[ average_surface_dir 'average_white_surface_2_left_',  num2str(NumMesh) '_t1.obj' ],  ...
        [ average_surface_dir 'average_white_surface_2_right_', num2str(NumMesh) '_t1.obj' ] });
    TemplateSurf{8} = SurfStatReadSurf({[ average_surface_dir 'average_white_surface_3_left_',  num2str(NumMesh) '_t1.obj' ],  ...
        [ average_surface_dir 'average_white_surface_3_right_', num2str(NumMesh) '_t1.obj' ] });
    TemplateSurf_standard = SurfStatReadSurf({'/local_raid/seokjun/01_project/04_IntracorticalAnalysis/01_analysis/surf_reg_model_left.obj', ...
                                              '/local_raid/seokjun/01_project/04_IntracorticalAnalysis/01_analysis/surf_reg_model_right.obj'});
    TemplateMask = SurfStatMaskCut(TemplateSurf_standard);
    AAL_data = SurfStatReadData1('/local_raid/seokjun/01_project/04_IntracorticalAnalysis/01_analysis/aal_both_rsl_final.txt');
   
end

%% Compute z-score
for zscore_computation = 1
    
    % Modality: rs-fMRI
    
    % Feature:
    %          ALFF_IntraCortical   (numControl x vertexnum),
    %          ReHo_IntraCortical   (numControl x vertexnum),
    %          DC_IntraCortical     (numControl x vertexnum),
    %          BC_IntraCortical     (numControl x vertexnum),
    
    % healthy controls
    for healthy_control_computation = 1
        
        % rsfMRI
        for rsfMRI_read_file = 1
            
            % Local functional features
            rsfMRI_ALFF_IntCortical_temp      = zeros(size(case_num_healthy_cont, 1), vertexnum);
            rsfMRI_ReHo_IntCortical_temp      = zeros(size(case_num_healthy_cont, 1), vertexnum);
            rsfMRI_DC_IntCortical_temp        = zeros(size(case_num_healthy_cont, 1), vertexnum);
            rsfMRI_BC_IntCortical_temp        = zeros(size(case_num_healthy_cont, 1), vertexnum);
            
            rsfMRI_ALFF_IntCortical_mean      = zeros(1, vertexnum);
            rsfMRI_ReHo_IntCortical_mean      = zeros(1, vertexnum);
            rsfMRI_DC_IntCortical_mean        = zeros(1, vertexnum);
            rsfMRI_BC_IntCortical_mean        = zeros(1, vertexnum);
            
            data_dir = '/host/gypsy/local_raid/seokjun/01_project/05_NoelRest/03_Result/01_DPARSF_processing/DPARSFA_scrub/02_without_mean_sig_regout/bbr_rest/';
            file_postfix = [ '_sm5_rsl.txt' ];
            for j = 1 : size(case_num_healthy_cont, 1)
                
                fprintf([ case_num_healthy_cont{j} ' : ' ]);
                prefix_path = [ data_dir '/' Prefix_healthy_cont '_' case_num_healthy_cont{j} '/' Prefix_healthy_cont '_' case_num_healthy_cont{j} ];
                
                %% ALFF. ReHo, DC, BC
                rsfMRI_ALFF_IntCortical_temp(j, 1:vertexnum) = SurfStatReadData( { [ prefix_path '_zALFFMap_mid_left_' num2str(NumMesh) file_postfix ], [ prefix_path '_zALFFMap_mid_right_' num2str(NumMesh) file_postfix ] } );
                rsfMRI_ReHo_IntCortical_temp(j, 1:vertexnum) = SurfStatReadData( { [ prefix_path '_zReHoMap_mid_left_' num2str(NumMesh) file_postfix ], [ prefix_path '_zReHoMap_mid_right_' num2str(NumMesh) file_postfix ] } );
                rsfMRI_DC_IntCortical_temp(j, 1:vertexnum)   = SurfStatReadData( { [ prefix_path '_zDegreeCentrality_PositiveBinarizedSumBrainMap_mid_left_' num2str(NumMesh) file_postfix ], [ prefix_path '_zDegreeCentrality_PositiveBinarizedSumBrainMap_mid_right_' num2str(NumMesh) file_postfix ] } );
                rsfMRI_BC_IntCortical_temp(j, 1:vertexnum)   = SurfStatReadData( { [ prefix_path '_zDegreeCentrality_PositiveBinarizedSumBrainMap_mid_left_' num2str(NumMesh) file_postfix ], [ prefix_path '_zDegreeCentrality_PositiveBinarizedSumBrainMap_mid_right_' num2str(NumMesh) file_postfix ] } );
                fprintf([ 'ALFF, ReHo, DC, BC in Mid\n' ]);
                
            end        
            
            %% Mean of functional features
            rsfMRI_ALFF_IntCortical_mean      = mean(rsfMRI_ALFF_IntCortical_temp(:, 1:vertexnum), 1);
            rsfMRI_ReHo_IntCortical_mean      = mean(rsfMRI_ReHo_IntCortical_temp(:, 1:vertexnum), 1);
            rsfMRI_DC_IntCortical_mean        = mean(rsfMRI_DC_IntCortical_temp(:, 1:vertexnum), 1);
            rsfMRI_BC_IntCortical_mean        = mean(rsfMRI_BC_IntCortical_temp(:, 1:vertexnum), 1);
            rsfMRI_ALFF_IntCortical_std       = std(rsfMRI_ALFF_IntCortical_temp(:, 1:vertexnum), 0, 1);
            rsfMRI_ReHo_IntCortical_std       = std(rsfMRI_ReHo_IntCortical_temp(:, 1:vertexnum), 0, 1);
            rsfMRI_DC_IntCortical_std         = std(rsfMRI_DC_IntCortical_temp(:, 1:vertexnum), 0, 1);
            rsfMRI_BC_IntCortical_std         = std(rsfMRI_BC_IntCortical_temp(:, 1:vertexnum), 0, 1);
            
            %% 1st mean and SD examination ...
            IntMask = repmat(TemplateMask, [5, 1]);
            SubMask = repmat(TemplateMask, [3, 1]);
            
        end
        
        % T1 Outlier detection
        % Total outlier: 306, 313, 322
        % compute mean and SD using rest of controls
        for compute_final_mean_std = 1
            
            outlier_cases = { '306_1', '313_1', '322_2' };
            case_cont_real = 1:size(case_num_healthy_cont, 1);
            outlier_vector = zeros(size(case_num_healthy_cont, 1), 1);
            for i = 1 : size(outlier_cases, 2)
                outlier_vector = outlier_vector + strcmp(case_num_healthy_cont, outlier_cases{i});
            end
            case_cont_real(logical(outlier_vector)) = [];
            
            %% Recalculate mean and SD of T1 features
            rsfMRI_ALFF_IntCortical_mean      = mean(rsfMRI_ALFF_IntCortical_temp(case_cont_real, 1:vertexnum), 1);
            rsfMRI_ReHo_IntCortical_mean      = mean(rsfMRI_ReHo_IntCortical_temp(case_cont_real, 1:vertexnum), 1);
            rsfMRI_DC_IntCortical_mean        = mean(rsfMRI_DC_IntCortical_temp(case_cont_real, 1:vertexnum), 1);
            rsfMRI_BC_IntCortical_mean        = mean(rsfMRI_BC_IntCortical_temp(case_cont_real, 1:vertexnum), 1);
            
            rsfMRI_ALFF_IntCortical_std       = std(rsfMRI_ALFF_IntCortical_temp(case_cont_real, 1:vertexnum), 0, 1);
            rsfMRI_ReHo_IntCortical_std       = std(rsfMRI_ReHo_IntCortical_temp(case_cont_real, 1:vertexnum), 0, 1);
            rsfMRI_DC_IntCortical_std         = std(rsfMRI_DC_IntCortical_temp(case_cont_real, 1:vertexnum), 0, 1);
            rsfMRI_BC_IntCortical_std         = std(rsfMRI_BC_IntCortical_temp(case_cont_real, 1:vertexnum), 0, 1);
        end
        
        for rsfMRI_zscore_healthy_cont = 1
            
            case_num_healthy_cont_temp = case_num_healthy_cont(case_cont_real); 
            
            %% z-score
            rsfMRI_ALFF_IntCortical_z_healthy_cont        = zeros(size(case_num_healthy_cont_temp, 1), vertexnum);
            rsfMRI_ReHo_IntCortical_z_healthy_cont        = zeros(size(case_num_healthy_cont_temp, 1), vertexnum);
            rsfMRI_DC_IntCortical_z_healthy_cont          = zeros(size(case_num_healthy_cont_temp, 1), vertexnum);
            rsfMRI_BC_IntCortical_z_healthy_cont          = zeros(size(case_num_healthy_cont_temp, 1), vertexnum);
            
            %% rsfMRI: calculate z-score of patients w.r.t control's distribution
            for j = 1 : size(case_num_healthy_cont_temp, 1)
                fprintf([ case_num_healthy_cont_temp{j} ' : ' ]);
                
                %% ALFF, ReHo, DC, BC
                rsfMRI_ALFF_IntCortical_z_healthy_cont(j, 1:vertexnum) = (rsfMRI_ALFF_IntCortical_temp(case_cont_real(j), 1:vertexnum) - rsfMRI_ALFF_IntCortical_mean) ./ rsfMRI_ALFF_IntCortical_std;
                rsfMRI_ReHo_IntCortical_z_healthy_cont(j, 1:vertexnum) = (rsfMRI_ReHo_IntCortical_temp(case_cont_real(j), 1:vertexnum) - rsfMRI_ReHo_IntCortical_mean) ./ rsfMRI_ReHo_IntCortical_std;
                rsfMRI_DC_IntCortical_z_healthy_cont(j, 1:vertexnum)   = (rsfMRI_DC_IntCortical_temp(case_cont_real(j), 1:vertexnum)   - rsfMRI_DC_IntCortical_mean)   ./ rsfMRI_DC_IntCortical_std;
                rsfMRI_BC_IntCortical_z_healthy_cont(j, 1:vertexnum)   = (rsfMRI_BC_IntCortical_temp(case_cont_real(j), 1:vertexnum)   - rsfMRI_BC_IntCortical_mean)   ./ rsfMRI_BC_IntCortical_std;                
                fprintf([ 'z-score of ALFF. ReHo, DC and BC in Mid\n' ]);
            end
            
        end
    end
    
    % disease controls
    for disease_control_computation = 1

        % rsfMRI
        for rsfMRI_read_file = 1
            
            % Local functional features
            rsfMRI_ALFF_IntCortical_temp2      = zeros(size(case_num_disease_cont, 1), vertexnum);
            rsfMRI_ReHo_IntCortical_temp2      = zeros(size(case_num_disease_cont, 1), vertexnum);
            rsfMRI_DC_IntCortical_temp2        = zeros(size(case_num_disease_cont, 1), vertexnum);
            rsfMRI_BC_IntCortical_temp2        = zeros(size(case_num_disease_cont, 1), vertexnum);            
            
            data_dir = '/host/weka/export02/data/min/fMRI/DPARSF_processing/DPARSFA_scrub/02_without_mean_sig_regout/bbr_rest/';
            file_postfix = [ '_sm5_rsl.txt' ];
            for j = 1 : size(case_num_disease_cont, 1)
                
                fprintf([ case_num_disease_cont{j} ' : ' ]);
                prefix_path = [ data_dir '/' Prefix_disease_cont '_' case_num_disease_cont{j} '/' Prefix_disease_cont '_' case_num_disease_cont{j} ];
                
                %% ALFF. ReHo, DC, BC
                rsfMRI_ALFF_IntCortical_temp2(j, 1:vertexnum) = SurfStatReadData( { [ prefix_path '_zALFFMap_mid_left_' num2str(NumMesh) file_postfix ], [ prefix_path '_zALFFMap_mid_right_' num2str(NumMesh) file_postfix ] } );
                rsfMRI_ReHo_IntCortical_temp2(j, 1:vertexnum) = SurfStatReadData( { [ prefix_path '_zReHoMap_mid_left_' num2str(NumMesh) file_postfix ], [ prefix_path '_zReHoMap_mid_right_' num2str(NumMesh) file_postfix ] } );
                rsfMRI_DC_IntCortical_temp2(j, 1:vertexnum)   = SurfStatReadData( { [ prefix_path '_zDegreeCentrality_PositiveBinarizedSumBrainMap_mid_left_' num2str(NumMesh) file_postfix ], [ prefix_path '_zDegreeCentrality_PositiveBinarizedSumBrainMap_mid_right_' num2str(NumMesh) file_postfix ] } );
                rsfMRI_BC_IntCortical_temp2(j, 1:vertexnum)   = SurfStatReadData( { [ prefix_path '_zDegreeCentrality_PositiveBinarizedSumBrainMap_mid_left_' num2str(NumMesh) file_postfix ], [ prefix_path '_zDegreeCentrality_PositiveBinarizedSumBrainMap_mid_right_' num2str(NumMesh) file_postfix ] } );
                fprintf([ 'ALFF, ReHo, DC, BC in Mid\n' ]);
                
            end        
            
        end

        for rsfMRI_zscore_disease_cont = 1
            
            %% z-score
            rsfMRI_ALFF_IntCortical_z_disease_cont        = zeros(size(case_num_disease_cont, 1), vertexnum);
            rsfMRI_ReHo_IntCortical_z_disease_cont        = zeros(size(case_num_disease_cont, 1), vertexnum);
            rsfMRI_DC_IntCortical_z_disease_cont          = zeros(size(case_num_disease_cont, 1), vertexnum);
            rsfMRI_BC_IntCortical_z_disease_cont          = zeros(size(case_num_disease_cont, 1), vertexnum);
            
            %% rsfMRI: calculate z-score of patients w.r.t control's distribution
            for j = 1 : size(case_num_disease_cont, 1)
                fprintf([ case_num_disease_cont{j} ' : ' ]);
                
                %% ALFF, ReHo, DC, BC
                rsfMRI_ALFF_IntCortical_z_disease_cont(j, 1:vertexnum) = (rsfMRI_ALFF_IntCortical_temp2(j, 1:vertexnum) - rsfMRI_ALFF_IntCortical_mean) ./ rsfMRI_ALFF_IntCortical_std;
                rsfMRI_ReHo_IntCortical_z_disease_cont(j, 1:vertexnum) = (rsfMRI_ReHo_IntCortical_temp2(j, 1:vertexnum) - rsfMRI_ReHo_IntCortical_mean) ./ rsfMRI_ReHo_IntCortical_std;
                rsfMRI_DC_IntCortical_z_disease_cont(j, 1:vertexnum)   = (rsfMRI_DC_IntCortical_temp2(j, 1:vertexnum) - rsfMRI_DC_IntCortical_mean) ./ rsfMRI_DC_IntCortical_std;
                rsfMRI_BC_IntCortical_z_disease_cont(j, 1:vertexnum)   = (rsfMRI_BC_IntCortical_temp2(j, 1:vertexnum) - rsfMRI_BC_IntCortical_mean) ./ rsfMRI_BC_IntCortical_std;                
                fprintf([ 'z-score of ALFF. ReHo, DC and BC in Mid\n' ]);
            end
            
        end

        for visualization_zscore = 1
            
            %% 1st mean and SD examination ...
            IntMask = repmat(TemplateMask, [5, 1]);
            SubMask = repmat(TemplateMask, [3, 1]);
            
            for j = 1 : size(case_num_disease_cont, 1)
                OUTPATH_temp = '/local_raid/seokjun/01_project/04_IntracorticalAnalysis/03_result/97_3T_TLE_HS/';
                idx = find(strcmp(case_num_disease_cont, case_num_disease_cont{j}));
                
                %     f = figure; CSFSurfStatView(rsfMRI_FA_IntCortical_z_disease_cont(:, 1:vertexnum, idx).*IntMask, TemplateSurf(1:5), 'INC', 'FA intra z-score');         CSFSurfStatViewColLim([-5 5]);
                %     exportfigbo(f,[OUTPATH_temp '/' case_num_disease_cont{j} '/' Prefix_disease_cont '_' case_num_disease_cont{j} '_z-score_rsfMRI_Intra_FA.png' ], 'png', 6); close(f);
                %     f = figure; CSFSurfStatView(rsfMRI_MD_IntCortical_z_disease_cont(:, 1:vertexnum, idx).*IntMask, TemplateSurf(1:5), 'INC', 'MD intra z-score');         CSFSurfStatViewColLim([-5 5]);
                %     exportfigbo(f,[OUTPATH_temp '/' case_num_disease_cont{j} '/' Prefix_disease_cont '_' case_num_disease_cont{j} '_z-score_rsfMRI_Intra_MD.png' ], 'png', 6); close(f);
                %     f = figure; CSFSurfStatView(rsfMRI_FA_SubCortical_z_disease_cont(1:5, 1:vertexnum, idx).*IntMask, TemplateSurf(1:5), 'INC', 'FA sub z-score');         CSFSurfStatViewColLim([-5 5]);
                %     exportfigbo(f,[OUTPATH_temp '/' case_num_disease_cont{j} '/' Prefix_disease_cont '_' case_num_disease_cont{j} '_z-score_rsfMRI_Sub_FA.png' ], 'png', 6); close(f);
                %     f = figure; CSFSurfStatView(rsfMRI_MD_SubCortical_z_disease_cont(1:5, 1:vertexnum, idx).*IntMask, TemplateSurf(1:5), 'INC', 'MD sub z-score');         CSFSurfStatViewColLim([-5 5]);
                %     exportfigbo(f,[OUTPATH_temp '/' case_num_disease_cont{j} '/' Prefix_disease_cont '_' case_num_disease_cont{j} '_z-score_rsfMRI_Sub_MD.png' ], 'png', 6); close(f);
                
            end
            
        end
        
    end
    
end

%% Extract lesion features and save as a mat file
for lesional_feature = 1
    
    % extract features from lesional areas using manually segmented labels
    % average the values from multiple vertices
    for lesion_feature_extraction = 1
        
        lesion_label_dir = '/local_raid/seokjun/CTFCD-1.2.0_64/Lesion/lesion_surf/';
        
        % disease controls        
        mean_z_lesion_disease_cont = zeros(4, size(case_num_pat, 1), size(case_num_disease_cont, 1));
        std_z_lesion_disease_cont  = zeros(4, size(case_num_pat, 1), size(case_num_disease_cont, 1));  
        case_num_pat_temp = case_num_pat;        
        for j = 1 : size(case_num_pat_temp, 1)
            idx = find(strcmp(case_num_pat_temp, case_num_pat_temp{j}));
            
            if(exist([ lesion_label_dir '/' Prefix_pat '_' case_num_pat_temp{j} '_label_union_left_rsl.txt' ], 'file'))
                lesion_label_data = SurfStatReadData( { [ lesion_label_dir '/' Prefix_pat '_' case_num_pat_temp{j} '_label_union_left_rsl.txt' ], ...
                                                        [ lesion_label_dir '/' Prefix_pat '_' case_num_pat_temp{j} '_label_union_right_rsl.txt' ] } );
                
                lesion_label_v = find(lesion_label_data ~= 0);
                
                for c = 1 : size(case_num_disease_cont, 1)
                    %% FA MD Intra/Subcortical surfaces
                    mean_z_lesion_disease_cont(1, j, c) = mean(rsfMRI_ALFF_IntCortical_z_disease_cont(c, lesion_label_v), 2);
                    std_z_lesion_disease_cont(1, j, c)  = std(rsfMRI_ALFF_IntCortical_z_disease_cont(c, lesion_label_v), 0, 2);
                    mean_z_lesion_disease_cont(2, j, c) = mean(rsfMRI_ReHo_IntCortical_z_disease_cont(c, lesion_label_v), 2);
                    std_z_lesion_disease_cont(2, j, c)  = std(rsfMRI_ReHo_IntCortical_z_disease_cont(c, lesion_label_v), 0, 2);
                    mean_z_lesion_disease_cont(3, j, c) = mean(rsfMRI_DC_IntCortical_z_disease_cont(c, lesion_label_v), 2);
                    std_z_lesion_disease_cont(3, j, c)  = std(rsfMRI_DC_IntCortical_z_disease_cont(c, lesion_label_v), 0, 2);
                    mean_z_lesion_disease_cont(4, j, c) = mean(rsfMRI_BC_IntCortical_z_disease_cont(c, lesion_label_v), 2);
                    std_z_lesion_disease_cont(4, j, c)  = std(rsfMRI_BC_IntCortical_z_disease_cont(c, lesion_label_v), 0, 2);
                end
            else
                disp([ case_num_pat_temp{j} ' does not have a lesion label file' ]);
            end
        end
        
        % healthy controls        
        mean_z_lesion_healthy_cont = zeros(4, size(case_num_pat, 1), size(case_num_healthy_cont_temp, 1));
        std_z_lesion_healthy_cont  = zeros(4, size(case_num_pat, 1), size(case_num_healthy_cont_temp, 1));  
        case_num_pat_temp = case_num_pat;        
        for j = 1 : size(case_num_pat_temp, 1)
            idx = find(strcmp(case_num_pat_temp, case_num_pat_temp{j}));
            
            if(exist([ lesion_label_dir '/' Prefix_pat '_' case_num_pat_temp{j} '_label_union_left_rsl.txt' ], 'file'))
                lesion_label_data = SurfStatReadData( { [ lesion_label_dir '/' Prefix_pat '_' case_num_pat_temp{j} '_label_union_left_rsl.txt' ], ...
                                                        [ lesion_label_dir '/' Prefix_pat '_' case_num_pat_temp{j} '_label_union_right_rsl.txt' ] } );
                
                lesion_label_v = find(lesion_label_data ~= 0);
                
                for c = 1 : size(case_num_healthy_cont_temp, 1)
                    %% FA MD Intra/Subcortical surfaces
                    mean_z_lesion_healthy_cont(1, j, c) = mean(rsfMRI_ALFF_IntCortical_z_healthy_cont(c, lesion_label_v), 2);
                    std_z_lesion_healthy_cont(1, j, c)  = std(rsfMRI_ALFF_IntCortical_z_healthy_cont(c, lesion_label_v), 0, 2);
                    mean_z_lesion_healthy_cont(2, j, c) = mean(rsfMRI_ReHo_IntCortical_z_healthy_cont(c, lesion_label_v), 2);
                    std_z_lesion_healthy_cont(2, j, c)  = std(rsfMRI_ReHo_IntCortical_z_healthy_cont(c, lesion_label_v), 0, 2);
                    mean_z_lesion_healthy_cont(3, j, c) = mean(rsfMRI_DC_IntCortical_z_healthy_cont(c, lesion_label_v), 2);
                    std_z_lesion_healthy_cont(3, j, c)  = std(rsfMRI_DC_IntCortical_z_healthy_cont(c, lesion_label_v), 0, 2);
                    mean_z_lesion_healthy_cont(4, j, c) = mean(rsfMRI_BC_IntCortical_z_healthy_cont(c, lesion_label_v), 2);
                    std_z_lesion_healthy_cont(4, j, c)  = std(rsfMRI_BC_IntCortical_z_healthy_cont(c, lesion_label_v), 0, 2);
                end
            else
                disp([ case_num_pat_temp{j} ' does not have a lesion label file' ]);
            end
        end
        
    end

    for compute_lesion_size = 1
        
        lesion_label_volume = [ 23385, 8276, 5400, 217, 332, 285, 1131, 6167, 1087, 3581, 143, 361, 14115, 5113, 203, 1417, 8119, 223, 847, 45, 585, 585, 832, 0, 2545, 1366, 8672, 417, 516, 264, 2679, 3305, 0, 0]';
        lesion_label_size = zeros(size(case_num_pat, 1), 1);
        for j = 1 : size(case_num_pat_temp, 1)
            idx = find(strcmp(case_num_pat_temp, case_num_pat_temp{j}));
            
            if(exist([ lesion_label_dir '/' Prefix_pat '_' case_num_pat_temp{j} '_label_union_left_rsl.txt' ], 'file'))
                lesion_label_data = SurfStatReadData( { [ lesion_label_dir '/' Prefix_pat '_' case_num_pat_temp{j} '_label_union_left_rsl.txt' ], ...
                    [ lesion_label_dir '/' Prefix_pat '_' case_num_pat_temp{j} '_label_union_right_rsl.txt' ] } );
                
                lesion_label_size(j) = sum(lesion_label_data ~= 0);
            end
        end
        
    end
    
    for save_zscore = 1
        
        if(~exist('zscore_database_rsfMRI_disease_cont.mat', 'file'))
            save('zscore_database_rsfMRI_disease_cont.mat', ... 
                 'rsfMRI_ALFF_IntCortical_z_disease_cont', 'rsfMRI_ReHo_IntCortical_z_disease_cont', 'rsfMRI_DC_IntCortical_z_disease_cont', 'rsfMRI_BC_IntCortical_z_disease_cont', ...
                 'mean_z_lesion_disease_cont', 'std_z_lesion_disease_cont');
        end
        
        if(~exist('zscore_database_rsfMRI_healthy_cont.mat', 'file'))
            save('zscore_database_rsfMRI_healthy_cont.mat', ... 
                 'rsfMRI_ALFF_IntCortical_z_healthy_cont', 'rsfMRI_ReHo_IntCortical_z_healthy_cont', 'rsfMRI_DC_IntCortical_z_healthy_cont', 'rsfMRI_BC_IntCortical_z_healthy_cont', ...
                 'mean_z_lesion_healthy_cont', 'std_z_lesion_healthy_cont');
        end
        
    end

end