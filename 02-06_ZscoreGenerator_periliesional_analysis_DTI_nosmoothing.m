clear all
close all
addpath('/data/noel/noel2/local/matlab/surfstat_chicago');

%% Read demograpic data
for read_demodata = 1
    
    OutDir                    = '/local_raid/seokjun/01_project/04_IntracorticalAnalysis/03_result/'
    Group_cont                = 'control'
    Prefix_cont               = 'TLE'
    Group_pat                 = 'FCD'
    Prefix_pat                = 'mcd'
    Cases_cont                = '/local_raid/seokjun/01_project/04_IntracorticalAnalysis/01_analysis/Demographic_data_control_DTI.txt'
    Cases_pat                 = '/local_raid/seokjun/01_project/04_IntracorticalAnalysis/01_analysis/Demographic_data_FCD_DTI.txt'
    
    Left_Right                = 'both'
    NumIntSurf                = 3
    NumSubSurf                = 5
    NumMesh                   = 81920
    SamplingSpace             = 'native'
    img_contrast              = {'t1'}
    Kernel                    = 2
    Parametric                = 'quadratic'
    average_surface_dir       = '/local_raid/seokjun/01_project/04_IntracorticalAnalysis/03_result/average_surfaces/'
    
    visualization = 1;
    load('/local_raid/seokjun/01_project/04_IntracorticalAnalysis/01_analysis/pipeline/postprocessing_zscore_visualization/colormap_noel1.mat');
    
    fid = fopen(Cases_cont);
    demo = textscan(fid, '%s%f%s', 'Delimiter', ',', 'CollectOutput', 1);
    case_num_cont = demo{1};
    age_cont = demo{2};
    gender_cont = demo{3};
    fclose(fid);
    
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
    
    mean_age_cont = mean(age_cont);
    std_age_cont = std(age_cont, 0);
    
    mean_age_pat = mean(age_pat);
    std_age_pat = std(age_pat, 0);
    
    [h,p,ci,stats] = ttest2(age_cont,age_pat);
    
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
    
    % Modality: DTI
    
    % Feature:
    %          FA_IntraCortical (vertexnum x 5 x numControl),
    %          FA_SubCortical   (vertexnum x 7 x numControl),
    %          MD_IntraCortical (vertexnum x 5 x numControl),
    %          MD_SubCortical   (vertexnum x 7 x numControl),
    
    % Abb: FA fractional anisotropy
    %      MD mean diffusivity
    % DTI index features: FA, MD
    
    % controls
    for control_computation = 1
        
        % DTI
        for DTI_read_file = 1
            
            DTI_FA_IntCortical_temp      = zeros(NumIntSurf+2, vertexnum, size(case_num_cont, 1));
            DTI_FA_SubCortical_temp      = zeros(NumSubSurf, vertexnum, size(case_num_cont, 1));
            DTI_MD_IntCortical_temp      = zeros(NumIntSurf+2, vertexnum, size(case_num_cont, 1));
            DTI_MD_SubCortical_temp      = zeros(NumSubSurf, vertexnum, size(case_num_cont, 1));
            
            DTI_FA_IntCortical_mean      = zeros(NumIntSurf+2, vertexnum);
            DTI_FA_SubCortical_mean      = zeros(NumSubSurf, vertexnum);
            DTI_MD_IntCortical_mean      = zeros(NumIntSurf+2, vertexnum);
            DTI_MD_SubCortical_mean      = zeros(NumSubSurf, vertexnum);
            DTI_FA_IntCortical_std       = zeros(NumIntSurf+2, vertexnum);
            DTI_FA_SubCortical_std       = zeros(NumSubSurf, vertexnum);
            DTI_MD_IntCortical_std       = zeros(NumIntSurf+2, vertexnum);
            DTI_MD_SubCortical_std       = zeros(NumSubSurf, vertexnum);
            
            data_dir = '/local_raid/seokjun/01_project/04_IntracorticalAnalysis/03_result/';
            file_postfix = [ '_rsl.txt' ];
            for j = 1 : size(case_num_cont, 1)
                fprintf([ case_num_cont{j} ' : ' ]);
                prefix_path = [ data_dir '/' case_num_cont{j} '/measurement/' Prefix_cont '_' case_num_cont{j} ];
                postfix_surf = '_native_dti_';
                
                %% Intensity-based features: corrected RI, pg, tg
                for i = 1 : NumIntSurf+2
                    if(i == 1)
                        basename{1}  = [ prefix_path '_gray_surface_left_'  num2str(NumMesh) postfix_surf ];
                        basename{2}  = [ prefix_path '_gray_surface_right_'  num2str(NumMesh) postfix_surf ];
                    elseif(i == 5)
                        basename{1}  = [ prefix_path '_white_surface_left_'  num2str(NumMesh) postfix_surf ];
                        basename{2}  = [ prefix_path '_white_surface_right_'  num2str(NumMesh) postfix_surf ];
                    else
                        basename{1}  = [ prefix_path '_intracortical_surface_' num2str(i-1) '_left_'  num2str(NumMesh) postfix_surf ];
                        basename{2}  = [ prefix_path '_intracortical_surface_' num2str(i-1) '_right_'  num2str(NumMesh) postfix_surf ];
                    end
                    
                    DTI_FA_IntCortical_temp(i, 1:vertexnum, j) = SurfStatReadData( { [ basename{1} 'FA' file_postfix ], [ basename{2} 'FA' file_postfix ] } );
                    DTI_MD_IntCortical_temp(i, 1:vertexnum, j) = SurfStatReadData( { [ basename{1} 'MD' file_postfix ], [ basename{2} 'MD' file_postfix ] } );
                    
                end
                fprintf([ 'FA, MD in Inc, ' ]);
                
                %% Intensity-based features of subcortical surfaces: RI, pg, tg
                for i = 1 : NumSubSurf
                    basename{1}  = [ prefix_path '_white_surface_' num2str(i) '_left_'  num2str(NumMesh) postfix_surf ];
                    basename{2}  = [ prefix_path '_white_surface_' num2str(i) '_right_'  num2str(NumMesh) postfix_surf ];
                    
                    DTI_FA_SubCortical_temp(i, 1:vertexnum, j) = SurfStatReadData( { [ basename{1} 'FA' file_postfix ], [ basename{2} 'FA' file_postfix ] } );
                    DTI_MD_SubCortical_temp(i, 1:vertexnum, j) = SurfStatReadData( { [ basename{1} 'MD' file_postfix ], [ basename{2} 'MD' file_postfix ] } );
                end
                fprintf([ 'FA, MD in Sub\n' ]);
            end
            
            %% Mean of T1 features
            DTI_FA_IntCortical_mean      = mean(DTI_FA_IntCortical_temp(:, 1:vertexnum, :), 3);
            DTI_FA_SubCortical_mean      = mean(DTI_FA_SubCortical_temp(:, 1:vertexnum, :), 3);
            DTI_MD_IntCortical_mean      = mean(DTI_MD_IntCortical_temp(:, 1:vertexnum, :), 3);
            DTI_MD_SubCortical_mean      = mean(DTI_MD_SubCortical_temp(:, 1:vertexnum, :), 3);
            DTI_FA_IntCortical_std       = std(DTI_FA_IntCortical_temp(:, 1:vertexnum, :), 0, 3);
            DTI_FA_SubCortical_std       = std(DTI_FA_SubCortical_temp(:, 1:vertexnum, :), 0, 3);
            DTI_MD_IntCortical_std       = std(DTI_MD_IntCortical_temp(:, 1:vertexnum, :), 0, 3);
            DTI_MD_SubCortical_std       = std(DTI_MD_SubCortical_temp(:, 1:vertexnum, :), 0, 3);
            
            %% 1st mean and SD examination ...
            IntMask = repmat(TemplateMask, [5, 1]);
            SubMask = repmat(TemplateMask, [3, 1]);
            
            % figure; CSFSurfStatView(DTI_FA_IntCortical_mean.*IntMask, TemplateSurf(1:5), 'INC', 'FA mean'); CSFSurfStatViewColLim([0 0.6]);
            % figure; CSFSurfStatView(DTI_FA_IntCortical_std.*IntMask,  TemplateSurf(1:5), 'INC', 'FA std');  CSFSurfStatViewColLim([0 0.3]);
            % figure; CSFSurfStatView(DTI_MD_IntCortical_mean.*IntMask, TemplateSurf(1:5), 'INC', 'MD mean'); CSFSurfStatViewColLim([0 0.002]);
            % figure; CSFSurfStatView(DTI_MD_IntCortical_std.*IntMask,  TemplateSurf(1:5), 'INC', 'MD std');  CSFSurfStatViewColLim([0 0.001]);
            % figure; CSFSurfStatView(DTI_FA_SubCortical_mean(1:5,:), TemplateSurf(1:5), 'INC', 'Sub FA mean'); CSFSurfStatViewColLim([0 0.8]);
            % figure; CSFSurfStatView(DTI_FA_SubCortical_std(1:5,:),  TemplateSurf(1:5), 'INC', 'Sub FA std');  CSFSurfStatViewColLim([0 0.3]);
            % figure; CSFSurfStatView(DTI_MD_SubCortical_mean(1:5,:), TemplateSurf(1:5), 'INC', 'Sub MD mean'); CSFSurfStatViewColLim([0.0005 0.001]);
            % figure; CSFSurfStatView(DTI_MD_SubCortical_std(1:5,:),  TemplateSurf(1:5), 'INC', 'Sub MD std');  CSFSurfStatViewColLim([0 0.0005]);
            
        end
        
        % Total outlier: 306, 313, 322
        % compute mean and SD using rest of controls
        for compute_final_mean_std = 1
            
            outlier_cases = { '306_1', '313_1', '322_1' };
            case_cont_real = 1:size(case_num_cont, 1);
            outlier_vector = zeros(size(case_num_cont, 1), 1);
            for i = 1 : size(outlier_cases, 2)
                outlier_vector = outlier_vector + strcmp(case_num_cont, outlier_cases{i});
            end
            case_cont_real(logical(outlier_vector)) = [];
            
            %% Recalculate mean and SD of T1 features
            DTI_FA_IntCortical_mean      = mean(DTI_FA_IntCortical_temp(:, 1:vertexnum, case_cont_real), 3);
            DTI_FA_SubCortical_mean      = mean(DTI_FA_SubCortical_temp(:, 1:vertexnum, case_cont_real), 3);
            DTI_MD_IntCortical_mean      = mean(DTI_MD_IntCortical_temp(:, 1:vertexnum, case_cont_real), 3);
            DTI_MD_SubCortical_mean      = mean(DTI_MD_SubCortical_temp(:, 1:vertexnum, case_cont_real), 3);
            
            DTI_FA_IntCortical_std       = std(DTI_FA_IntCortical_temp(:, 1:vertexnum, case_cont_real), 0, 3);
            DTI_FA_SubCortical_std       = std(DTI_FA_SubCortical_temp(:, 1:vertexnum, case_cont_real), 0, 3);
            DTI_MD_IntCortical_std       = std(DTI_MD_IntCortical_temp(:, 1:vertexnum, case_cont_real), 0, 3);
            DTI_MD_SubCortical_std       = std(DTI_MD_SubCortical_temp(:, 1:vertexnum, case_cont_real), 0, 3);
            
            
        end
        
    end
    
    
    % patients
    for patient_computation = 1
        
        % read DTI files
        for DTI_read_file = 1
            
            %% DTI index features: FA, MD
            
            DTI_FA_IntCortical_temp2      = zeros(NumIntSurf+2, vertexnum, size(case_num_cont, 1));
            DTI_FA_SubCortical_temp2      = zeros(NumSubSurf, vertexnum, size(case_num_cont, 1));
            DTI_MD_IntCortical_temp2      = zeros(NumIntSurf+2, vertexnum, size(case_num_cont, 1));
            DTI_MD_SubCortical_temp2      = zeros(NumSubSurf, vertexnum, size(case_num_cont, 1));
            
            DTI_FA_IntCortical_mean2      = zeros(NumIntSurf+2, vertexnum);
            DTI_FA_SubCortical_mean2      = zeros(NumSubSurf, vertexnum);
            DTI_MD_IntCortical_mean2      = zeros(NumIntSurf+2, vertexnum);
            DTI_MD_SubCortical_mean2      = zeros(NumSubSurf, vertexnum);
            DTI_FA_IntCortical_std2       = zeros(NumIntSurf+2, vertexnum);
            DTI_FA_SubCortical_std2       = zeros(NumSubSurf, vertexnum);
            DTI_MD_IntCortical_std2       = zeros(NumIntSurf+2, vertexnum);
            DTI_MD_SubCortical_std2       = zeros(NumSubSurf, vertexnum);
            
            data_dir = '/local_raid/seokjun/01_project/04_IntracorticalAnalysis/03_result/';
            file_postfix = [ '_rsl.txt' ];
            for j = 1 : size(case_num_pat, 1)
                fprintf([ case_num_pat{j} ' : ' ]);
                
                prefix_path = [ data_dir '/' case_num_pat{j} '/measurement/' Prefix_pat '_' case_num_pat{j} ];
                postfix_surf = '_native_dti_';
                
                %% Intensity-based features: corrected RI, pg, tg
                for i = 1 : NumIntSurf+2
                    if(i == 1)
                        basename{1}  = [ prefix_path '_gray_surface_left_'  num2str(NumMesh) postfix_surf ];
                        basename{2}  = [ prefix_path '_gray_surface_right_'  num2str(NumMesh) postfix_surf ];
                    elseif(i == 5)
                        basename{1}  = [ prefix_path '_white_surface_left_'  num2str(NumMesh) postfix_surf ];
                        basename{2}  = [ prefix_path '_white_surface_right_'  num2str(NumMesh) postfix_surf ];
                    else
                        basename{1}  = [ prefix_path '_intracortical_surface_' num2str(i-1) '_left_'  num2str(NumMesh) postfix_surf ];
                        basename{2}  = [ prefix_path '_intracortical_surface_' num2str(i-1) '_right_'  num2str(NumMesh) postfix_surf ];
                    end
                    
                    DTI_FA_IntCortical_temp2(i, 1:vertexnum, j) = SurfStatReadData( { [ basename{1} 'FA' file_postfix ], [ basename{2} 'FA' file_postfix ] } );
                    DTI_MD_IntCortical_temp2(i, 1:vertexnum, j) = SurfStatReadData( { [ basename{1} 'MD' file_postfix ], [ basename{2} 'MD' file_postfix ] } );
                    
                end
                fprintf([ 'FA, MD in Inc, ' ]);
                
                %% Intensity-based features of subcortical surfaces: RI, pg, tg
                for i = 1 : NumSubSurf
                    basename{1}  = [ prefix_path '_white_surface_' num2str(i) '_left_'  num2str(NumMesh) postfix_surf ];
                    basename{2}  = [ prefix_path '_white_surface_' num2str(i) '_right_'  num2str(NumMesh) postfix_surf ];
                    
                    DTI_FA_SubCortical_temp2(i, 1:vertexnum, j) = SurfStatReadData( { [ basename{1} 'FA' file_postfix ], [ basename{2} 'FA' file_postfix ] } );
                    DTI_MD_SubCortical_temp2(i, 1:vertexnum, j) = SurfStatReadData( { [ basename{1} 'MD' file_postfix ], [ basename{2} 'MD' file_postfix ] } );
                end
                fprintf([ 'FA, MD in Sub\n' ]);
                
            end
            
        end
        
        % compute DTI z-score
        for DTI_zscore = 1
            
            %% z-score
            DTI_FA_IntCortical_z         = zeros(NumIntSurf+2, vertexnum, size(case_num_pat, 1));
            DTI_FA_SubCortical_z         = zeros(NumSubSurf,   vertexnum, size(case_num_pat, 1));
            DTI_MD_IntCortical_z         = zeros(NumIntSurf+2, vertexnum, size(case_num_pat, 1));
            DTI_MD_SubCortical_z         = zeros(NumSubSurf,   vertexnum, size(case_num_pat, 1));
            
            %% DTI: calculate z-score of patients w.r.t control's distribution
            for j = 1 : size(case_num_pat, 1)
                fprintf([ case_num_pat{j} ' : ' ]);
                
                %% Intensity-based features: corrected RI, pg, tg
                DTI_FA_IntCortical_z(:, 1:vertexnum, j) = (DTI_FA_IntCortical_temp2(:, 1:vertexnum, j) - DTI_FA_IntCortical_mean) ./ DTI_FA_IntCortical_std;
                DTI_MD_IntCortical_z(:, 1:vertexnum, j) = (DTI_MD_IntCortical_temp2(:, 1:vertexnum, j) - DTI_MD_IntCortical_mean) ./ DTI_MD_IntCortical_std;
                fprintf([ 'z-score of FA and MD in Inc, ' ]);
                
                %% Intensity-based features of subcortical surfaces: RI, pg, tg
                DTI_FA_SubCortical_z(:, 1:vertexnum, j) = (DTI_FA_SubCortical_temp2(:, 1:vertexnum, j) - DTI_FA_SubCortical_mean) ./ DTI_FA_SubCortical_std;
                DTI_MD_SubCortical_z(:, 1:vertexnum, j) = (DTI_MD_SubCortical_temp2(:, 1:vertexnum, j) - DTI_MD_SubCortical_mean) ./ DTI_MD_SubCortical_std;
                fprintf([ 'z-score of FA and MD in Sub\n' ]);
            end
            
        end
        
        for visualization_zscore = 1
            
            %% 1st mean and SD examination ...
            IntMask = repmat(TemplateMask, [5, 1]);
            SubMask = repmat(TemplateMask, [3, 1]);
            
            for j = 1 : size(case_num_pat, 1)
                OUTPATH_temp = '/local_raid/seokjun/01_project/04_IntracorticalAnalysis/03_result/';
                idx = find(strcmp(case_num_pat, case_num_pat{j}));
                
                f = figure; CSFSurfStatView(DTI_FA_IntCortical_z(:, 1:vertexnum, idx).*IntMask, TemplateSurf(1:5), 'INC', 'FA intra z-score');         CSFSurfStatViewColLim([-5 5]);
                exportfigbo(f,[OUTPATH_temp '/' case_num_pat{j} '/' Prefix_pat '_' case_num_pat{j} '_z-score_DTI_Intra_FA.png' ], 'png', 6); close(f);
                f = figure; CSFSurfStatView(DTI_MD_IntCortical_z(:, 1:vertexnum, idx).*IntMask, TemplateSurf(1:5), 'INC', 'MD intra z-score');         CSFSurfStatViewColLim([-5 5]);
                exportfigbo(f,[OUTPATH_temp '/' case_num_pat{j} '/' Prefix_pat '_' case_num_pat{j} '_z-score_DTI_Intra_MD.png' ], 'png', 6); close(f);
                f = figure; CSFSurfStatView(DTI_FA_SubCortical_z(1:5, 1:vertexnum, idx).*IntMask, TemplateSurf(1:5), 'INC', 'FA sub z-score');         CSFSurfStatViewColLim([-5 5]);
                exportfigbo(f,[OUTPATH_temp '/' case_num_pat{j} '/' Prefix_pat '_' case_num_pat{j} '_z-score_DTI_Sub_FA.png' ], 'png', 6); close(f);
                f = figure; CSFSurfStatView(DTI_MD_SubCortical_z(1:5, 1:vertexnum, idx).*IntMask, TemplateSurf(1:5), 'INC', 'MD sub z-score');         CSFSurfStatViewColLim([-5 5]);
                exportfigbo(f,[OUTPATH_temp '/' case_num_pat{j} '/' Prefix_pat '_' case_num_pat{j} '_z-score_DTI_Sub_MD.png' ], 'png', 6); close(f);
            end
            
        end
        
    end
    
end
