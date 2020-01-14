clear all
close all
addpath('/data/noel/noel2/local/matlab/surfstat_chicago');

%% Read demograpic data
for read_demodata = 1
    
    OutDir                    = '/local_raid/seokjun/01_project/04_IntracorticalAnalysis/03_result/'
    
    Cases_pat                 = '/local_raid/seokjun/01_project/04_IntracorticalAnalysis/01_analysis/Demographic_data_FCD.txt'
    Group_pat                 = 'FCD'
    Prefix_pat                = 'mcd'
    
    Cases_healthy_cont        = '/local_raid/seokjun/01_project/04_IntracorticalAnalysis/01_analysis/Demographic_data_control.txt'    
    Group_healthy_cont        = 'control'
    Prefix_healthy_cont       = 'TLE'
    
    Cases_disease_cont        = '/local_raid/seokjun/01_project/04_IntracorticalAnalysis/01_analysis/Demographic_data_TLE_final.txt'
    Group_disease_cont        = 'TLE'
    Prefix_disease_cont       = 'TLE'
        
    Left_Right                = 'both'
    NumIntSurf                = 3
    NumSubSurf                = 3
    NumMesh                   = 81920
    SamplingSpace             = 'native'
    img_contrast              = {'t1'}
    Kernel                    = 5
    Parametric                = 'quadratic'
    average_surface_dir       = '/local_raid/seokjun/01_project/04_IntracorticalAnalysis/03_result/average_surfaces/'
    
    visualization = 1;
    load('/local_raid/seokjun/01_project/04_IntracorticalAnalysis/01_analysis/pipeline/postprocessing_zscore_visualization/colormap_noel1.mat');
    
    fid                         = fopen(Cases_pat);
    demo                        = textscan(fid, '%s%f%s%d%s%d%d%s%s', 'Delimiter', ',', 'CollectOutput', 1);
    case_num_pat                = demo{1};
    age_pat                     = demo{2};
    gender_pat                  = demo{3};
    lesion_volume               = demo{4};
    histo_type                  = demo{5};
    initial                     = demo{6}(:, 1);
    transmantle                 = demo{6}(:, 2);
    location                    = demo{7}(:, 1);
    seizure_lateralization      = demo{7}(:, 2);
    fclose(fid);
    
    fid                         = fopen(Cases_healthy_cont);
    demo                        = textscan(fid, '%s%f%s', 'Delimiter', ',', 'CollectOutput', 1);
    case_num_healthy_cont       = demo{1};
    age_healthy_cont            = demo{2};
    gender_healthy_cont         = demo{3};
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
    
    % Modality: T1, FLAIR, IR
    
    % Feature: RI_IntraCortical (vertexnum x 5 x numControl),
    %          RI_SubCortical   (vertexnum x 3 x numControl),
    %          PG_IntraCortical (vertexnum x 5 x numControl),
    %          PG_GM_WM         (vertexnum x numControl),
    %          PG_SubCortical   (vertexnum x 3 x numControl),
    %          TG_IntraCortical (vertexnum x 5 x numControl),
    %          TG_SubCortical   (vertexnum x 3 x numControl),
    %          CT, MC, SD       (vertexnum x numControl, only T1)
    
    % Abb: RI relative intensity
    %      PG perpendicular gradient
    %      TG tangential gradient
    %      CT cortical thickness
    %      MC mean curvature
    %      SD sulcal depth
    
    % healthy controls
    for healthy_control_computation = 1
        
        % T1
        for T1_read_file = 1
            
            % Intensity-based features
            T1_RI_IntCortical_temp      = zeros(NumIntSurf+2, vertexnum, size(case_num_healthy_cont, 1));
            T1_RI_SubCortical_temp      = zeros(NumSubSurf, vertexnum, size(case_num_healthy_cont, 1));
            T1_PG_IntCortical_temp      = zeros(NumIntSurf+2, vertexnum, size(case_num_healthy_cont, 1));
            T1_PG_gw_IntCortical_temp   = zeros(size(case_num_healthy_cont, 1), vertexnum);
            T1_PG_SubCortical_temp      = zeros(NumSubSurf, vertexnum, size(case_num_healthy_cont, 1));
            T1_TG_IntCortical_temp      = zeros(NumIntSurf+2, vertexnum, size(case_num_healthy_cont, 1));
            T1_TG_SubCortical_temp      = zeros(NumSubSurf, vertexnum, size(case_num_healthy_cont, 1));
            
            % Morphological features
            T1_CT_midCortical_temp      = zeros(size(case_num_healthy_cont, 1), vertexnum);
            T1_SD_wmCortical_temp       = zeros(size(case_num_healthy_cont, 1), vertexnum);
            T1_MC_midCortical_temp      = zeros(size(case_num_healthy_cont, 1), vertexnum);
            
            T1_RI_IntCortical_mean      = zeros(NumIntSurf+2, vertexnum);
            T1_RI_SubCortical_mean      = zeros(NumSubSurf, vertexnum);
            T1_PG_IntCortical_mean      = zeros(NumIntSurf+2, vertexnum);
            T1_PG_SubCortical_mean      = zeros(NumSubSurf, vertexnum);
            T1_TG_IntCortical_mean      = zeros(NumIntSurf+2, vertexnum);
            T1_TG_SubCortical_mean      = zeros(NumSubSurf, vertexnum);
            T1_PG_gw_IntCortical_mean   = zeros(1, vertexnum);
            T1_CT_midCortical_mean      = zeros(1, vertexnum);
            T1_MC_midCortical_mean      = zeros(1, vertexnum);
            T1_SD_wmCortical_mean       = zeros(1, vertexnum);
            T1_RI_IntCortical_std       = zeros(NumIntSurf+2, vertexnum);
            T1_RI_SubCortical_std       = zeros(NumSubSurf, vertexnum);
            T1_PG_IntCortical_std       = zeros(NumIntSurf+2, vertexnum);
            T1_PG_SubCortical_std       = zeros(NumSubSurf, vertexnum);
            T1_TG_IntCortical_std       = zeros(NumIntSurf+2, vertexnum);
            T1_TG_SubCortical_std       = zeros(NumSubSurf, vertexnum);
            T1_PG_gw_IntCortical_std    = zeros(1, vertexnum);
            T1_CT_midCortical_std       = zeros(1, vertexnum);
            T1_MC_midCortical_std       = zeros(1, vertexnum);
            T1_SD_wmCortical_std        = zeros(1, vertexnum);
            
            data_dir = '/local_raid/seokjun/01_project/04_IntracorticalAnalysis/03_result/';
            file_postfix = [ '_quadratic_sm_5_rsl.txt' ];
            for j = 1 : size(case_num_healthy_cont, 1)
                fprintf([ case_num_healthy_cont{j} ' : ' ]);
                prefix_path = [ data_dir '/' case_num_healthy_cont{j} '/measurement/' Prefix_healthy_cont '_' case_num_healthy_cont{j} ];
                postfix_surf = '_native_t1_';
                
                % Intensity-based features: corrected RI, pg, tg
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
                    
                    T1_RI_IntCortical_temp(i, 1:vertexnum, j) = SurfStatReadData( { [ basename{1} 'RI_corrected' file_postfix ], [ basename{2} 'RI_corrected' file_postfix ] } );
                    T1_PG_IntCortical_temp(i, 1:vertexnum, j) = SurfStatReadData( { [ basename{1} 'pg' file_postfix ], [ basename{2} 'pg' file_postfix ] } );
                    T1_TG_IntCortical_temp(i, 1:vertexnum, j) = SurfStatReadData( { [ basename{1} 'tg' file_postfix ], [ basename{2} 'tg' file_postfix ] } );
                end
                fprintf([ 'T1: RI_corrected, PG, TG in Inc, ' ]);
                
                % Morphological / Intensity-based features: MC, SD, PG_GW
                basename{1}  = [ prefix_path '_white_surface_left_'  num2str(NumMesh) postfix_surf ];
                basename{2}  = [ prefix_path '_white_surface_right_'  num2str(NumMesh) postfix_surf ];
                T1_PG_gw_IntCortical_temp(j, 1:vertexnum) = SurfStatReadData( { [ basename{1} 'pg_GM_WM' file_postfix ], [ basename{2} 'pg_GM_WM' file_postfix ] } );
                
                basename{1}  = [ prefix_path '_intracortical_surface_2_left_'  num2str(NumMesh) postfix_surf ];
                basename{2}  = [ prefix_path '_intracortical_surface_2_right_'  num2str(NumMesh) postfix_surf ];
                T1_CT_midCortical_temp(j, 1:vertexnum) = SurfStatReadData( { [ basename{1} 'ct' file_postfix ], [ basename{2} 'ct' file_postfix ] } );
                
                basename{1}  = [ prefix_path '_intracortical_surface_2_left_'  num2str(NumMesh) postfix_surf ];
                basename{2}  = [ prefix_path '_intracortical_surface_2_right_'  num2str(NumMesh) postfix_surf ];
                T1_MC_midCortical_temp(j, 1:vertexnum) = SurfStatReadData( { [ basename{1} 'mc' file_postfix ], [ basename{2} 'mc' file_postfix ] } );
                
                basename{1}  = [ prefix_path '_white_surface_left_'  num2str(NumMesh) postfix_surf ];
                basename{2}  = [ prefix_path '_white_surface_right_'  num2str(NumMesh) postfix_surf ];
                T1_SD_wmCortical_temp(j, 1:vertexnum) = SurfStatReadData( { [ basename{1} 'sd2' file_postfix ], [ basename{2} 'sd2' file_postfix ] } );
                fprintf([ 'CT, MC, SD, PG_gw in Inc, ' ]);
                
                % Intensity-based features of subcortical surfaces: RI, pg, tg
                for i = 1 : NumSubSurf
                    basename{1}  = [ prefix_path '_white_surface_' num2str(i) '_left_'  num2str(NumMesh) postfix_surf ];
                    basename{2}  = [ prefix_path '_white_surface_' num2str(i) '_right_'  num2str(NumMesh) postfix_surf ];
                    
                    T1_RI_SubCortical_temp(i, 1:vertexnum, j) = SurfStatReadData( { [ basename{1} 'RI' file_postfix ], [ basename{2} 'RI' file_postfix ] } );
                    T1_PG_SubCortical_temp(i, 1:vertexnum, j) = SurfStatReadData( { [ basename{1} 'pg' file_postfix ], [ basename{2} 'pg' file_postfix ] } );
                    T1_TG_SubCortical_temp(i, 1:vertexnum, j) = SurfStatReadData( { [ basename{1} 'tg' file_postfix ], [ basename{2} 'tg' file_postfix ] } );
                end
                fprintf([ 'RI, PG, TG in Sub\n' ]);
            end
            
            %% Mean of T1 features
            T1_RI_IntCortical_mean      = mean(T1_RI_IntCortical_temp(:, 1:vertexnum, :), 3);
            T1_RI_SubCortical_mean      = mean(T1_RI_SubCortical_temp(:, 1:vertexnum, :), 3);
            T1_PG_IntCortical_mean      = mean(T1_PG_IntCortical_temp(:, 1:vertexnum, :), 3);
            T1_PG_SubCortical_mean      = mean(T1_PG_SubCortical_temp(:, 1:vertexnum, :), 3);
            T1_TG_IntCortical_mean      = mean(T1_TG_IntCortical_temp(:, 1:vertexnum, :), 3);
            T1_TG_SubCortical_mean      = mean(T1_TG_SubCortical_temp(:, 1:vertexnum, :), 3);
            T1_PG_gw_IntCortical_mean   = mean(T1_PG_gw_IntCortical_temp(:, 1:vertexnum), 1);
            T1_CT_midCortical_mean      = mean(T1_CT_midCortical_temp(:, 1:vertexnum), 1);
            T1_MC_midCortical_mean      = mean(T1_MC_midCortical_temp(:, 1:vertexnum), 1);
            T1_SD_wmCortical_mean       = mean(T1_SD_wmCortical_temp(:, 1:vertexnum), 1);
            
            T1_RI_IntCortical_std       = std(T1_RI_IntCortical_temp(:, 1:vertexnum, :), 0, 3);
            T1_RI_SubCortical_std       = std(T1_RI_SubCortical_temp(:, 1:vertexnum, :), 0, 3);
            T1_PG_IntCortical_std       = std(T1_PG_IntCortical_temp(:, 1:vertexnum, :), 0, 3);
            T1_PG_SubCortical_std       = std(T1_PG_SubCortical_temp(:, 1:vertexnum, :), 0, 3);
            T1_TG_IntCortical_std       = std(T1_TG_IntCortical_temp(:, 1:vertexnum, :), 0, 3);
            T1_TG_SubCortical_std       = std(T1_TG_SubCortical_temp(:, 1:vertexnum, :), 0, 3);
            T1_PG_gw_IntCortical_std    = std(T1_PG_gw_IntCortical_temp(:, 1:vertexnum), 0, 1);
            T1_CT_midCortical_std       = std(T1_CT_midCortical_temp(:, 1:vertexnum), 0, 1);
            T1_MC_midCortical_std       = std(T1_MC_midCortical_temp(:, 1:vertexnum), 0, 1);
            T1_SD_wmCortical_std        = std(T1_SD_wmCortical_temp(:, 1:vertexnum), 0, 1);
            
            %% 1st mean and SD examination ...
            IntMask = repmat(TemplateMask, [5, 1]);
            SubMask = repmat(TemplateMask, [3, 1]);
            
            % figure; CSFSurfStatView(T1_RI_IntCortical_mean.*IntMask, TemplateSurf(1:5), 'INC', 'RI mean'); CSFSurfStatViewColLim([-200 200]);
            % figure; CSFSurfStatView(T1_RI_IntCortical_std.*IntMask,  TemplateSurf(1:5), 'INC', 'RI std');  CSFSurfStatViewColLim([0 50]);
            % figure; CSFSurfStatView(T1_PG_IntCortical_mean.*IntMask, TemplateSurf(1:5), 'INC', 'PG mean'); CSFSurfStatViewColLim([0 300]);
            % figure; CSFSurfStatView(T1_PG_IntCortical_std.*IntMask,  TemplateSurf(1:5), 'INC', 'PG std');  CSFSurfStatViewColLim([0 150]);
            % figure; CSFSurfStatView(T1_TG_IntCortical_mean.*IntMask, TemplateSurf(1:5), 'INC', 'TG mean'); CSFSurfStatViewColLim([0 50]);
            % figure; CSFSurfStatView(T1_TG_IntCortical_std.*IntMask,  TemplateSurf(1:5), 'INC', 'TG std');  CSFSurfStatViewColLim([0 20]);
            % figure; CSFSurfStatView(T1_PG_gw_IntCortical_mean.*TemplateMask, TemplateSurf{3}, 'MF', 'PG mean'); CSFSurfStatViewColLim([0 300]);
            % figure; CSFSurfStatView(T1_PG_gw_IntCortical_std.*TemplateMask,  TemplateSurf{3}, 'MF', 'PG std');  CSFSurfStatViewColLim([0 150]);
            % figure; CSFSurfStatView(T1_CT_midCortical_mean.*TemplateMask, TemplateSurf{3}, 'MF', 'CT mean'); CSFSurfStatViewColLim([0 5]);
            % figure; CSFSurfStatView(T1_CT_midCortical_std.*TemplateMask,  TemplateSurf{3}, 'MF', 'CT std');  CSFSurfStatViewColLim([0 1]);
            % figure; CSFSurfStatView(T1_MC_midCortical_mean.*TemplateMask, TemplateSurf{3}, 'MF', 'MC mean'); CSFSurfStatViewColLim([0 0.3]);
            % figure; CSFSurfStatView(T1_MC_midCortical_std.*TemplateMask,  TemplateSurf{3}, 'MF', 'MC std');  CSFSurfStatViewColLim([0 0.15]);
            % figure; CSFSurfStatView(T1_SD_wmCortical_mean.*TemplateMask, TemplateSurf{3},  'MF', 'SD mean'); CSFSurfStatViewColLim([0 25]);
            % figure; CSFSurfStatView(T1_SD_wmCortical_std.*TemplateMask,  TemplateSurf{3},  'MF', 'SD std');  CSFSurfStatViewColLim([0 10]);
            %
            % figure; CSFSurfStatView(T1_RI_SubCortical_mean.*SubMask, TemplateSurf(6:8), 'SBC', 'RI mean'); CSFSurfStatViewColLim([-100 100]);
            % figure; CSFSurfStatView(T1_RI_SubCortical_std.*SubMask,  TemplateSurf(6:8), 'SBC', 'RI std');  CSFSurfStatViewColLim([0 100]);
            % figure; CSFSurfStatView(T1_PG_SubCortical_mean.*SubMask, TemplateSurf(6:8), 'SBC', 'PG mean'); CSFSurfStatViewColLim([0 100]);
            % figure; CSFSurfStatView(T1_PG_SubCortical_std.*SubMask,  TemplateSurf(6:8), 'SBC', 'PG std');  CSFSurfStatViewColLim([0 50]);
            % figure; CSFSurfStatView(T1_TG_SubCortical_mean.*SubMask, TemplateSurf(6:8), 'SBC', 'TG mean'); CSFSurfStatViewColLim([0 40]);
            % figure; CSFSurfStatView(T1_TG_SubCortical_std.*SubMask,  TemplateSurf(6:8), 'SBC', 'TG std');  CSFSurfStatViewColLim([0 15]);
            
        end
        
        % T1 Outlier detection
        for T1_outlier_detection = 1
            sid = 1;
            temp_mean = T1_RI_IntCortical_mean(sid, :);
            temp_std  = T1_RI_IntCortical_std(sid, :);
            temp_mean_new = temp_mean;
            temp_std_new = temp_std;
            
            vid_max = 100;
            sub_idx = zeros(size(case_num_healthy_cont, 1), 28);
            
            %% RI outliers
            for j = 1 : 5
                sid = j;
                for i = 1 : vid_max
                    disp([ num2str(j) ' surf, ' num2str(i) ' vertex' ]);
                    vid = i;
                    [b,idx,outliers] = deleteoutliers(T1_RI_IntCortical_temp(sid, vid, :), 0.1);
                    sub_idx(idx, j) = sub_idx(idx, j) + 1;
                end
            end
            
            %% PG outliers
            for j = 6 : 10
                sid = j-5;
                for i = 1 : vid_max
                    disp([ num2str(sid) ' surf, ' num2str(i) ' vertex' ]);
                    vid = i;
                    [b,idx,outliers] = deleteoutliers(T1_PG_IntCortical_temp(sid, vid, :), 0.1);
                    sub_idx(idx, j) = sub_idx(idx, j) + 1;
                end
            end
            
            %% TG outliers
            for j = 11 : 15
                sid = j-10;
                for i = 1 : vid_max
                    disp([ num2str(sid) ' surf, ' num2str(i) ' vertex' ]);
                    vid = i;
                    [b,idx,outliers] = deleteoutliers(T1_TG_IntCortical_temp(sid, vid, :), 0.1);
                    sub_idx(idx, j) = sub_idx(idx, j) + 1;
                end
            end
            
            %% CT,MC,SD,PG_gw outliers
            j = 16;
            for i = 1 : vid_max
                disp([ num2str(j) ' surf, ' num2str(i) ' vertex' ]);
                vid = i;
                [b,idx,outliers] = deleteoutliers(T1_CT_midCortical_temp(:, vid), 0.1);
                sub_idx(idx, j) = sub_idx(idx, j) + 1;
            end
            j = 17;
            for i = 1 : vid_max
                disp([ num2str(j) ' surf, ' num2str(i) ' vertex' ]);
                vid = i;
                [b,idx,outliers] = deleteoutliers(T1_PG_gw_IntCortical_temp(:, vid), 0.1);
                sub_idx(idx, j) = sub_idx(idx, j) + 1;
            end
            j = 18;
            for i = 1 : vid_max
                disp([ num2str(j) ' surf, ' num2str(i) ' vertex' ]);
                vid = i;
                [b,idx,outliers] = deleteoutliers(T1_MC_midCortical_temp(:, vid), 0.1);
                sub_idx(idx, j) = sub_idx(idx, j) + 1;
            end
            j = 19;
            for i = 1 : vid_max
                disp([ num2str(j) ' surf, ' num2str(i) ' vertex' ]);
                vid = i;
                [b,idx,outliers] = deleteoutliers(T1_SD_wmCortical_temp(:, vid), 0.1);
                sub_idx(idx, j) = sub_idx(idx, j) + 1;
            end
            
            %% RI outliers
            for j = 20 : 22
                sid = j-19;
                for i = 1 : vid_max
                    disp([ num2str(sid) ' surf, ' num2str(i) ' vertex' ]);
                    vid = i;
                    [b,idx,outliers] = deleteoutliers(T1_RI_SubCortical_temp(sid, vid, :), 0.1);
                    sub_idx(idx, j) = sub_idx(idx, j) + 1;
                end
            end
            
            %% PG outliers
            for j = 23 : 25
                sid = j-22;
                for i = 1 : vid_max
                    disp([ num2str(sid) ' surf, ' num2str(i) ' vertex' ]);
                    vid = i;
                    [b,idx,outliers] = deleteoutliers(T1_PG_SubCortical_temp(sid, vid, :), 0.1);
                    sub_idx(idx, j) = sub_idx(idx, j) + 1;
                end
            end
            
            %% TG outliers
            for j = 26 : 28
                sid = j-25;
                for i = 1 : vid_max
                    disp([ num2str(sid) ' surf, ' num2str(i) ' vertex' ]);
                    vid = i;
                    [b,idx,outliers] = deleteoutliers(T1_TG_SubCortical_temp(sid, vid, :), 0.1);
                    sub_idx(idx, j) = sub_idx(idx, j) + 1;
                end
            end
            
            [outlier_score, idx_outlier ] = sort(sum(sub_idx, 2)/2800, 'descend');
            strvcat(case_num_healthy_cont{idx_outlier})
            
            % T1 outlier: 322, 306
            
        end
        
        % FLAIR
        for FLAIR_read_file = 1
            
            % Intensity-based features
            FLAIR_RI_IntCortical_temp      = zeros(NumIntSurf+2, vertexnum, size(case_num_healthy_cont, 1));
            FLAIR_RI_SubCortical_temp      = zeros(NumSubSurf, vertexnum, size(case_num_healthy_cont, 1));
            FLAIR_PG_IntCortical_temp      = zeros(NumIntSurf+2, vertexnum, size(case_num_healthy_cont, 1));
            FLAIR_PG_gw_IntCortical_temp   = zeros(size(case_num_healthy_cont, 1), vertexnum);
            FLAIR_PG_SubCortical_temp      = zeros(NumSubSurf, vertexnum, size(case_num_healthy_cont, 1));
            FLAIR_TG_IntCortical_temp      = zeros(NumIntSurf+2, vertexnum, size(case_num_healthy_cont, 1));
            FLAIR_TG_SubCortical_temp      = zeros(NumSubSurf, vertexnum, size(case_num_healthy_cont, 1));
            
            FLAIR_RI_IntCortical_mean      = zeros(NumIntSurf+2, vertexnum);
            FLAIR_RI_SubCortical_mean      = zeros(NumSubSurf, vertexnum);
            FLAIR_PG_IntCortical_mean      = zeros(NumIntSurf+2, vertexnum);
            FLAIR_PG_SubCortical_mean      = zeros(NumSubSurf, vertexnum);
            FLAIR_TG_IntCortical_mean      = zeros(NumIntSurf+2, vertexnum);
            FLAIR_TG_SubCortical_mean      = zeros(NumSubSurf, vertexnum);
            FLAIR_PG_gw_IntCortical_mean   = zeros(1, vertexnum);
            FLAIR_RI_IntCortical_std       = zeros(NumIntSurf+2, vertexnum);
            FLAIR_RI_SubCortical_std       = zeros(NumSubSurf, vertexnum);
            FLAIR_PG_IntCortical_std       = zeros(NumIntSurf+2, vertexnum);
            FLAIR_PG_SubCortical_std       = zeros(NumSubSurf, vertexnum);
            FLAIR_TG_IntCortical_std       = zeros(NumIntSurf+2, vertexnum);
            FLAIR_TG_SubCortical_std       = zeros(NumSubSurf, vertexnum);
            FLAIR_PG_gw_IntCortical_std    = zeros(1, vertexnum);
            
            data_dir = '/local_raid/seokjun/01_project/04_IntracorticalAnalysis/03_result/';
            file_postfix = [ '_quadratic_sm_5_rsl.txt' ];
            for j = 1 : size(case_num_healthy_cont, 1)
                fprintf([ case_num_healthy_cont{j} ' : ' ]);
                prefix_path = [ data_dir '/' case_num_healthy_cont{j} '/measurement/' Prefix_healthy_cont '_' case_num_healthy_cont{j} ];
                postfix_surf = '_native_flair_';
                
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
                    
                    FLAIR_RI_IntCortical_temp(i, 1:vertexnum, j) = SurfStatReadData( { [ basename{1} 'RI_corrected' file_postfix ], [ basename{2} 'RI_corrected' file_postfix ] } );
                    FLAIR_PG_IntCortical_temp(i, 1:vertexnum, j) = SurfStatReadData( { [ basename{1} 'pg' file_postfix ], [ basename{2} 'pg' file_postfix ] } );
                    FLAIR_TG_IntCortical_temp(i, 1:vertexnum, j) = SurfStatReadData( { [ basename{1} 'tg' file_postfix ], [ basename{2} 'tg' file_postfix ] } );
                end
                
                basename{1}  = [ prefix_path '_white_surface_left_'  num2str(NumMesh) postfix_surf ];
                basename{2}  = [ prefix_path '_white_surface_right_'  num2str(NumMesh) postfix_surf ];
                FLAIR_PG_gw_IntCortical_temp(j, 1:vertexnum) = SurfStatReadData( { [ basename{1} 'pg_GM_WM' file_postfix ], [ basename{2} 'pg_GM_WM' file_postfix ] } );
                
                fprintf([ 'FLAIR: RI_corrected, PG, TG in Inc, ' ]);
                
                %% Intensity-based features of subcortical surfaces: RI, pg, tg
                for i = 1 : NumSubSurf
                    basename{1}  = [ prefix_path '_white_surface_' num2str(i) '_left_'  num2str(NumMesh) postfix_surf ];
                    basename{2}  = [ prefix_path '_white_surface_' num2str(i) '_right_'  num2str(NumMesh) postfix_surf ];
                    
                    FLAIR_RI_SubCortical_temp(i, 1:vertexnum, j) = SurfStatReadData( { [ basename{1} 'RI' file_postfix ], [ basename{2} 'RI' file_postfix ] } );
                    FLAIR_PG_SubCortical_temp(i, 1:vertexnum, j) = SurfStatReadData( { [ basename{1} 'pg' file_postfix ], [ basename{2} 'pg' file_postfix ] } );
                    FLAIR_TG_SubCortical_temp(i, 1:vertexnum, j) = SurfStatReadData( { [ basename{1} 'tg' file_postfix ], [ basename{2} 'tg' file_postfix ] } );
                end
                fprintf([ 'RI, PG, TG in Sub\n' ]);
            end
            
            %% Mean of FLAIR features
            FLAIR_RI_IntCortical_mean      = mean(FLAIR_RI_IntCortical_temp(:, 1:vertexnum, :), 3);
            FLAIR_RI_SubCortical_mean      = mean(FLAIR_RI_SubCortical_temp(:, 1:vertexnum, :), 3);
            FLAIR_PG_IntCortical_mean      = mean(FLAIR_PG_IntCortical_temp(:, 1:vertexnum, :), 3);
            FLAIR_PG_SubCortical_mean      = mean(FLAIR_PG_SubCortical_temp(:, 1:vertexnum, :), 3);
            FLAIR_TG_IntCortical_mean      = mean(FLAIR_TG_IntCortical_temp(:, 1:vertexnum, :), 3);
            FLAIR_TG_SubCortical_mean      = mean(FLAIR_TG_SubCortical_temp(:, 1:vertexnum, :), 3);
            FLAIR_PG_gw_IntCortical_mean   = mean(FLAIR_PG_gw_IntCortical_temp(:, 1:vertexnum), 1);
            
            FLAIR_RI_IntCortical_std       = std(FLAIR_RI_IntCortical_temp(:, 1:vertexnum, :), 0, 3);
            FLAIR_RI_SubCortical_std       = std(FLAIR_RI_SubCortical_temp(:, 1:vertexnum, :), 0, 3);
            FLAIR_PG_IntCortical_std       = std(FLAIR_PG_IntCortical_temp(:, 1:vertexnum, :), 0, 3);
            FLAIR_PG_SubCortical_std       = std(FLAIR_PG_SubCortical_temp(:, 1:vertexnum, :), 0, 3);
            FLAIR_TG_IntCortical_std       = std(FLAIR_TG_IntCortical_temp(:, 1:vertexnum, :), 0, 3);
            FLAIR_TG_SubCortical_std       = std(FLAIR_TG_SubCortical_temp(:, 1:vertexnum, :), 0, 3);
            FLAIR_PG_gw_IntCortical_std    = std(FLAIR_PG_gw_IntCortical_temp(:, 1:vertexnum), 0, 1);
            
            %% 1st mean and SD examination ...
            % figure; CSFSurfStatView(FLAIR_RI_IntCortical_mean.*IntMask, TemplateSurf(1:5), 'INC', 'RI mean'); CSFSurfStatViewColLim([-500 100]);
            % figure; CSFSurfStatView(FLAIR_RI_IntCortical_std.*IntMask,  TemplateSurf(1:5), 'INC', 'RI std');  CSFSurfStatViewColLim([0 200]);
            % figure; CSFSurfStatView(FLAIR_PG_IntCortical_mean.*IntMask, TemplateSurf(1:5), 'INC', 'PG mean'); CSFSurfStatViewColLim([0 600]);
            % figure; CSFSurfStatView(FLAIR_PG_IntCortical_std.*IntMask,  TemplateSurf(1:5), 'INC', 'PG std');  CSFSurfStatViewColLim([0 300]);
            % figure; CSFSurfStatView(FLAIR_TG_IntCortical_mean.*IntMask, TemplateSurf(1:5), 'INC', 'TG mean'); CSFSurfStatViewColLim([0 50]);
            % figure; CSFSurfStatView(FLAIR_TG_IntCortical_std.*IntMask,  TemplateSurf(1:5), 'INC', 'TG std');  CSFSurfStatViewColLim([0 20]);
            % figure; CSFSurfStatView(FLAIR_PG_gw_IntCortical_mean.*TemplateMask, TemplateSurf{3}, 'MF', 'PG mean'); CSFSurfStatViewColLim([0 300]);
            % figure; CSFSurfStatView(FLAIR_PG_gw_IntCortical_std.*TemplateMask,  TemplateSurf{3}, 'MF', 'PG std');  CSFSurfStatViewColLim([0 150]);
            %
            % figure; CSFSurfStatView(FLAIR_RI_SubCortical_mean.*SubMask, TemplateSurf(6:8), 'SBC', 'RI mean'); CSFSurfStatViewColLim([-200 200]);
            % figure; CSFSurfStatView(FLAIR_RI_SubCortical_std.*SubMask,  TemplateSurf(6:8), 'SBC', 'RI std');  CSFSurfStatViewColLim([0 200]);
            % figure; CSFSurfStatView(FLAIR_PG_SubCortical_mean.*SubMask, TemplateSurf(6:8), 'SBC', 'PG mean'); CSFSurfStatViewColLim([0 100]);
            % figure; CSFSurfStatView(FLAIR_PG_SubCortical_std.*SubMask,  TemplateSurf(6:8), 'SBC', 'PG std');  CSFSurfStatViewColLim([0 50]);
            % figure; CSFSurfStatView(FLAIR_TG_SubCortical_mean.*SubMask, TemplateSurf(6:8), 'SBC', 'TG mean'); CSFSurfStatViewColLim([0 40]);
            % figure; CSFSurfStatView(FLAIR_TG_SubCortical_std.*SubMask,  TemplateSurf(6:8), 'SBC', 'TG std');  CSFSurfStatViewColLim([0 15]);
            
        end
        
        % FLAIR Outlier detection
        for FLAIR_outlier_detection = 1
            
            sid = 1;
            temp_mean = FLAIR_RI_IntCortical_mean(sid, :);
            temp_std  = FLAIR_RI_IntCortical_std(sid, :);
            temp_mean_new = temp_mean;
            temp_std_new = temp_std;
            
            vid_max = 100;
            sub_idx = zeros(size(case_num_healthy_cont, 1), 25);
            
            %% RI outliers
            for j = 1 : 5
                sid = j;
                for i = 1 : vid_max
                    disp([ num2str(j) ' surf, ' num2str(i) ' vertex' ]);
                    vid = i;
                    [b,idx,outliers] = deleteoutliers(FLAIR_RI_IntCortical_temp(sid, vid, :), 0.1);
                    sub_idx(idx, j) = sub_idx(idx, j) + 1;
                end
            end
            
            %% PG outliers
            for j = 6 : 10
                sid = j-5;
                for i = 1 : vid_max
                    disp([ num2str(sid) ' surf, ' num2str(i) ' vertex' ]);
                    vid = i;
                    [b,idx,outliers] = deleteoutliers(FLAIR_PG_IntCortical_temp(sid, vid, :), 0.1);
                    sub_idx(idx, j) = sub_idx(idx, j) + 1;
                end
            end
            
            %% TG outliers
            for j = 11 : 15
                sid = j-10;
                for i = 1 : vid_max
                    disp([ num2str(sid) ' surf, ' num2str(i) ' vertex' ]);
                    vid = i;
                    [b,idx,outliers] = deleteoutliers(FLAIR_TG_IntCortical_temp(sid, vid, :), 0.1);
                    sub_idx(idx, j) = sub_idx(idx, j) + 1;
                end
            end
            
            %% PG_gw outliers
            j = 16;
            for i = 1 : vid_max
                disp([ num2str(j) ' surf, ' num2str(i) ' vertex' ]);
                vid = i;
                [b,idx,outliers] = deleteoutliers(FLAIR_PG_gw_IntCortical_temp(:, vid), 0.1);
                sub_idx(idx, j) = sub_idx(idx, j) + 1;
            end
            
            %% RI outliers
            for j = 17 : 19
                sid = j-16;
                for i = 1 : vid_max
                    disp([ num2str(sid) ' surf, ' num2str(i) ' vertex' ]);
                    vid = i;
                    [b,idx,outliers] = deleteoutliers(FLAIR_RI_SubCortical_temp(sid, vid, :), 0.1);
                    sub_idx(idx, j) = sub_idx(idx, j) + 1;
                end
            end
            
            %% PG outliers
            for j = 20 : 22
                sid = j-19;
                for i = 1 : vid_max
                    disp([ num2str(sid) ' surf, ' num2str(i) ' vertex' ]);
                    vid = i;
                    [b,idx,outliers] = deleteoutliers(FLAIR_PG_SubCortical_temp(sid, vid, :), 0.1);
                    sub_idx(idx, j) = sub_idx(idx, j) + 1;
                end
            end
            
            %% TG outliers
            for j = 23 : 25
                sid = j-22;
                for i = 1 : vid_max
                    disp([ num2str(sid) ' surf, ' num2str(i) ' vertex' ]);
                    vid = i;
                    [b,idx,outliers] = deleteoutliers(FLAIR_TG_SubCortical_temp(sid, vid, :), 0.1);
                    sub_idx(idx, j) = sub_idx(idx, j) + 1;
                end
            end
            
            [outlier_score, idx_outlier ] = sort(sum(sub_idx, 2)/2800, 'descend');
            case_num_healthy_cont{idx_outlier}
            
            %% FLAIR outlier: 313, 306
            
        end
        
        % Total outlier: 306, 313, 322
        % compute mean and SD using rest of controls
        for compute_final_mean_std = 1
            
            outlier_cases = { '306_1', '313_1', '322_1' };
            case_cont_real = 1:size(case_num_healthy_cont, 1);
            outlier_vector = zeros(size(case_num_healthy_cont, 1), 1);
            for i = 1 : size(outlier_cases, 2)
                outlier_vector = outlier_vector + strcmp(case_num_healthy_cont, outlier_cases{i});
            end
            case_cont_real(logical(outlier_vector)) = [];
            
            %% Recalculate mean and SD of T1 features
            T1_RI_IntCortical_mean      = mean(T1_RI_IntCortical_temp(:, 1:vertexnum, case_cont_real), 3);
            T1_RI_SubCortical_mean      = mean(T1_RI_SubCortical_temp(:, 1:vertexnum, case_cont_real), 3);
            T1_PG_IntCortical_mean      = mean(T1_PG_IntCortical_temp(:, 1:vertexnum, case_cont_real), 3);
            T1_PG_SubCortical_mean      = mean(T1_PG_SubCortical_temp(:, 1:vertexnum, case_cont_real), 3);
            T1_TG_IntCortical_mean      = mean(T1_TG_IntCortical_temp(:, 1:vertexnum, case_cont_real), 3);
            T1_TG_SubCortical_mean      = mean(T1_TG_SubCortical_temp(:, 1:vertexnum, case_cont_real), 3);
            T1_PG_gw_IntCortical_mean   = mean(T1_PG_gw_IntCortical_temp(case_cont_real, 1:vertexnum), 1);
            T1_CT_midCortical_mean      = mean(T1_CT_midCortical_temp(case_cont_real, 1:vertexnum), 1);
            T1_MC_midCortical_mean      = mean(T1_MC_midCortical_temp(case_cont_real, 1:vertexnum), 1);
            T1_SD_wmCortical_mean       = mean(T1_SD_wmCortical_temp(case_cont_real, 1:vertexnum), 1);
            
            T1_RI_IntCortical_std       = std(T1_RI_IntCortical_temp(:, 1:vertexnum, case_cont_real), 0, 3);
            T1_RI_SubCortical_std       = std(T1_RI_SubCortical_temp(:, 1:vertexnum, case_cont_real), 0, 3);
            T1_PG_IntCortical_std       = std(T1_PG_IntCortical_temp(:, 1:vertexnum, case_cont_real), 0, 3);
            T1_PG_SubCortical_std       = std(T1_PG_SubCortical_temp(:, 1:vertexnum, case_cont_real), 0, 3);
            T1_TG_IntCortical_std       = std(T1_TG_IntCortical_temp(:, 1:vertexnum, case_cont_real), 0, 3);
            T1_TG_SubCortical_std       = std(T1_TG_SubCortical_temp(:, 1:vertexnum, case_cont_real), 0, 3);
            T1_PG_gw_IntCortical_std    = std(T1_PG_gw_IntCortical_temp(case_cont_real, 1:vertexnum), 0, 1);
            T1_CT_midCortical_std       = std(T1_CT_midCortical_temp(case_cont_real, 1:vertexnum), 0, 1);
            T1_MC_midCortical_std       = std(T1_MC_midCortical_temp(case_cont_real, 1:vertexnum), 0, 1);
            T1_SD_wmCortical_std        = std(T1_SD_wmCortical_temp(case_cont_real, 1:vertexnum), 0, 1);
            
            %% And FLAIR features too,
            FLAIR_RI_IntCortical_mean      = mean(FLAIR_RI_IntCortical_temp(:, 1:vertexnum, case_cont_real), 3);
            FLAIR_RI_SubCortical_mean      = mean(FLAIR_RI_SubCortical_temp(:, 1:vertexnum, case_cont_real), 3);
            FLAIR_PG_IntCortical_mean      = mean(FLAIR_PG_IntCortical_temp(:, 1:vertexnum, case_cont_real), 3);
            FLAIR_PG_SubCortical_mean      = mean(FLAIR_PG_SubCortical_temp(:, 1:vertexnum, case_cont_real), 3);
            FLAIR_TG_IntCortical_mean      = mean(FLAIR_TG_IntCortical_temp(:, 1:vertexnum, case_cont_real), 3);
            FLAIR_TG_SubCortical_mean      = mean(FLAIR_TG_SubCortical_temp(:, 1:vertexnum, case_cont_real), 3);
            FLAIR_PG_gw_IntCortical_mean   = mean(FLAIR_PG_gw_IntCortical_temp(case_cont_real, 1:vertexnum), 1);
            
            FLAIR_RI_IntCortical_std       = std(FLAIR_RI_IntCortical_temp(:, 1:vertexnum, case_cont_real), 0, 3);
            FLAIR_RI_SubCortical_std       = std(FLAIR_RI_SubCortical_temp(:, 1:vertexnum, case_cont_real), 0, 3);
            FLAIR_PG_IntCortical_std       = std(FLAIR_PG_IntCortical_temp(:, 1:vertexnum, case_cont_real), 0, 3);
            FLAIR_PG_SubCortical_std       = std(FLAIR_PG_SubCortical_temp(:, 1:vertexnum, case_cont_real), 0, 3);
            FLAIR_TG_IntCortical_std       = std(FLAIR_TG_IntCortical_temp(:, 1:vertexnum, case_cont_real), 0, 3);
            FLAIR_TG_SubCortical_std       = std(FLAIR_TG_SubCortical_temp(:, 1:vertexnum, case_cont_real), 0, 3);
            FLAIR_PG_gw_IntCortical_std    = std(FLAIR_PG_gw_IntCortical_temp(case_cont_real, 1:vertexnum), 0, 1);
            
        end
        
        % compute T1 z-score
        for T1_zscore_healthy_cont = 1
            
            case_num_cont_temp = case_num_healthy_cont(case_cont_real);            
            
            %% z-score
            T1_RI_IntCortical_z_healthy_cont         = zeros(NumIntSurf+2, vertexnum, size(case_num_cont_temp, 1));
            T1_RI_SubCortical_z_healthy_cont         = zeros(NumSubSurf,   vertexnum, size(case_num_cont_temp, 1));
            T1_PG_IntCortical_z_healthy_cont         = zeros(NumIntSurf+2, vertexnum, size(case_num_cont_temp, 1));
            T1_PG_SubCortical_z_healthy_cont         = zeros(NumSubSurf,   vertexnum, size(case_num_cont_temp, 1));
            T1_TG_IntCortical_z_healthy_cont         = zeros(NumIntSurf+2, vertexnum, size(case_num_cont_temp, 1));
            T1_TG_SubCortical_z_healthy_cont         = zeros(NumSubSurf,   vertexnum, size(case_num_cont_temp, 1));
            T1_CT_midCortical_z_healthy_cont         = zeros(size(case_num_cont_temp, 1), vertexnum);
            T1_MC_midCortical_z_healthy_cont         = zeros(size(case_num_cont_temp, 1), vertexnum);
            T1_SD_wmCortical_z_healthy_cont          = zeros(size(case_num_cont_temp, 1), vertexnum);
            T1_PG_gw_IntCortical_z_healthy_cont      = zeros(size(case_num_cont_temp, 1), vertexnum);
            
            %% T1: calculate z-score of patients w.r.t control's distribution
            for j = 1 : size(case_num_cont_temp, 1)
                fprintf([ case_num_cont_temp{j} ' : ' ]);
                
                %% Intensity-based features: corrected RI, pg, tg
                T1_RI_IntCortical_z_healthy_cont(:, 1:vertexnum, j) = (T1_RI_IntCortical_temp(:, 1:vertexnum, case_cont_real(j)) - T1_RI_IntCortical_mean) ./ T1_RI_IntCortical_std;
                T1_PG_IntCortical_z_healthy_cont(:, 1:vertexnum, j) = (T1_PG_IntCortical_temp(:, 1:vertexnum, case_cont_real(j)) - T1_PG_IntCortical_mean) ./ T1_PG_IntCortical_std;
                T1_TG_IntCortical_z_healthy_cont(:, 1:vertexnum, j) = (T1_TG_IntCortical_temp(:, 1:vertexnum, case_cont_real(j)) - T1_TG_IntCortical_mean) ./ T1_TG_IntCortical_std;
                T1_RI_IntCortical_z_healthy_cont(isnan(T1_RI_IntCortical_z_healthy_cont)) = 0; T1_RI_IntCortical_z_healthy_cont(isinf(T1_RI_IntCortical_z_healthy_cont)) = 0;
                T1_PG_IntCortical_z_healthy_cont(isnan(T1_PG_IntCortical_z_healthy_cont)) = 0; T1_PG_IntCortical_z_healthy_cont(isinf(T1_PG_IntCortical_z_healthy_cont)) = 0;
                T1_TG_IntCortical_z_healthy_cont(isnan(T1_TG_IntCortical_z_healthy_cont)) = 0; T1_PG_IntCortical_z_healthy_cont(isinf(T1_PG_IntCortical_z_healthy_cont)) = 0;
                fprintf([ 'z-score of T1: RI_corrected, PG, TG in Inc, ' ]);
                
                %% Morphological / Intensity-based features: MC, SD, PG_GW
                T1_PG_gw_IntCortical_z_healthy_cont(j, 1:vertexnum) = (T1_PG_gw_IntCortical_temp(case_cont_real(j), 1:vertexnum) - T1_PG_gw_IntCortical_mean) ./ T1_PG_gw_IntCortical_std;
                T1_CT_midCortical_z_healthy_cont(j, 1:vertexnum)    = (T1_CT_midCortical_temp(case_cont_real(j), 1:vertexnum) - T1_CT_midCortical_mean) ./ T1_CT_midCortical_std;
                T1_MC_midCortical_z_healthy_cont(j, 1:vertexnum)    = (T1_MC_midCortical_temp(case_cont_real(j), 1:vertexnum) - T1_MC_midCortical_mean) ./ T1_MC_midCortical_std;
                T1_SD_wmCortical_z_healthy_cont(j, 1:vertexnum)     = (T1_SD_wmCortical_temp(case_cont_real(j), 1:vertexnum) - T1_SD_wmCortical_mean) ./ T1_SD_wmCortical_std;
                T1_PG_gw_IntCortical_z_healthy_cont(isnan(T1_PG_gw_IntCortical_z_healthy_cont)) = 0; T1_PG_gw_IntCortical_z_healthy_cont(isinf(T1_PG_gw_IntCortical_z_healthy_cont)) = 0;
                T1_CT_midCortical_z_healthy_cont(isnan(T1_CT_midCortical_z_healthy_cont)) = 0; T1_CT_midCortical_z_healthy_cont(isinf(T1_CT_midCortical_z_healthy_cont)) = 0;
                T1_MC_midCortical_z_healthy_cont(isnan(T1_MC_midCortical_z_healthy_cont)) = 0; T1_MC_midCortical_z_healthy_cont(isinf(T1_MC_midCortical_z_healthy_cont)) = 0;
                T1_SD_wmCortical_z_healthy_cont(isnan(T1_SD_wmCortical_z_healthy_cont)) = 0; T1_SD_wmCortical_z_healthy_cont(isinf(T1_SD_wmCortical_z_healthy_cont)) = 0;
                fprintf([ 'CT, MC, SD, PG_gw in Inc, ' ]);
                
                %% Intensity-based features of subcortical surfaces: RI, pg, tg
                T1_RI_SubCortical_z_healthy_cont(:, 1:vertexnum, j) = (T1_RI_SubCortical_temp(:, 1:vertexnum, case_cont_real(j)) - T1_RI_SubCortical_mean) ./ T1_RI_SubCortical_std;
                T1_PG_SubCortical_z_healthy_cont(:, 1:vertexnum, j) = (T1_PG_SubCortical_temp(:, 1:vertexnum, case_cont_real(j)) - T1_PG_SubCortical_mean) ./ T1_PG_SubCortical_std;
                T1_TG_SubCortical_z_healthy_cont(:, 1:vertexnum, j) = (T1_TG_SubCortical_temp(:, 1:vertexnum, case_cont_real(j)) - T1_TG_SubCortical_mean) ./ T1_TG_SubCortical_std;
                T1_RI_SubCortical_z_healthy_cont(isnan(T1_RI_SubCortical_z_healthy_cont)) = 0; T1_RI_SubCortical_z_healthy_cont(isinf(T1_RI_SubCortical_z_healthy_cont)) = 0;
                T1_PG_SubCortical_z_healthy_cont(isnan(T1_PG_SubCortical_z_healthy_cont)) = 0; T1_PG_SubCortical_z_healthy_cont(isinf(T1_PG_SubCortical_z_healthy_cont)) = 0;
                T1_TG_SubCortical_z_healthy_cont(isnan(T1_TG_SubCortical_z_healthy_cont)) = 0; T1_TG_SubCortical_z_healthy_cont(isinf(T1_TG_SubCortical_z_healthy_cont)) = 0;
                fprintf([ 'RI, PG, TG in Sub\n' ]);
            end
            
        end

        % compute FLAIR z-score
        for FLAIR_zscore_healthy_cont = 1
            
            FLAIR_RI_IntCortical_z_healthy_cont         = zeros(NumIntSurf+2, vertexnum, size(case_num_cont_temp, 1));
            FLAIR_RI_SubCortical_z_healthy_cont         = zeros(NumSubSurf,   vertexnum, size(case_num_cont_temp, 1));
            FLAIR_PG_IntCortical_z_healthy_cont         = zeros(NumIntSurf+2, vertexnum, size(case_num_cont_temp, 1));
            FLAIR_PG_SubCortical_z_healthy_cont         = zeros(NumSubSurf,   vertexnum, size(case_num_cont_temp, 1));
            FLAIR_TG_IntCortical_z_healthy_cont         = zeros(NumIntSurf+2, vertexnum, size(case_num_cont_temp, 1));
            FLAIR_TG_SubCortical_z_healthy_cont         = zeros(NumSubSurf,   vertexnum, size(case_num_cont_temp, 1));
            FLAIR_PG_gw_IntCortical_z_healthy_cont      = zeros(size(case_num_disease_cont, 1), vertexnum);
            
            %% FLAIR: calculate z-score of patients w.r.t control's distribution
            for j = 1 : size(case_num_cont_temp, 1)
                fprintf([ case_num_cont_temp{j} ' : ' ]);
                
                %% Intensity-based features: corrected RI, pg, tg
                FLAIR_RI_IntCortical_z_healthy_cont(:, 1:vertexnum, j) = (FLAIR_RI_IntCortical_temp(:, 1:vertexnum, case_cont_real(j)) - FLAIR_RI_IntCortical_mean) ./ FLAIR_RI_IntCortical_std;
                FLAIR_PG_IntCortical_z_healthy_cont(:, 1:vertexnum, j) = (FLAIR_PG_IntCortical_temp(:, 1:vertexnum, case_cont_real(j)) - FLAIR_PG_IntCortical_mean) ./ FLAIR_PG_IntCortical_std;
                FLAIR_TG_IntCortical_z_healthy_cont(:, 1:vertexnum, j) = (FLAIR_TG_IntCortical_temp(:, 1:vertexnum, case_cont_real(j)) - FLAIR_TG_IntCortical_mean) ./ FLAIR_TG_IntCortical_std;
                FLAIR_RI_IntCortical_z_healthy_cont(isnan(FLAIR_RI_IntCortical_z_healthy_cont)) = 0; FLAIR_RI_IntCortical_z_healthy_cont(isinf(FLAIR_RI_IntCortical_z_healthy_cont)) = 0;
                FLAIR_PG_IntCortical_z_healthy_cont(isnan(FLAIR_PG_IntCortical_z_healthy_cont)) = 0; FLAIR_PG_IntCortical_z_healthy_cont(isinf(FLAIR_PG_IntCortical_z_healthy_cont)) = 0;
                FLAIR_TG_IntCortical_z_healthy_cont(isnan(FLAIR_TG_IntCortical_z_healthy_cont)) = 0; FLAIR_TG_IntCortical_z_healthy_cont(isinf(FLAIR_TG_IntCortical_z_healthy_cont)) = 0;
                fprintf([ 'z-score of FLAIR: RI_corrected, PG, TG in Inc, ' ]);
                
                %% Morphological / Intensity-based features: MC, SD, PG_GW
                FLAIR_PG_gw_IntCortical_z_healthy_cont(j, 1:vertexnum) = (FLAIR_PG_gw_IntCortical_temp(case_cont_real(j), 1:vertexnum) - FLAIR_PG_gw_IntCortical_mean) ./ FLAIR_PG_gw_IntCortical_std;
                FLAIR_PG_gw_IntCortical_z_healthy_cont(isnan(FLAIR_PG_gw_IntCortical_z_healthy_cont)) = 0;
                fprintf([ 'PG_gw in Inc, ' ]);
                
                %% Intensity-based features of subcortical surfaces: RI, pg, tg
                FLAIR_RI_SubCortical_z_healthy_cont(:, 1:vertexnum, j) = (FLAIR_RI_SubCortical_temp(:, 1:vertexnum, case_cont_real(j)) - FLAIR_RI_SubCortical_mean) ./ FLAIR_RI_SubCortical_std;
                FLAIR_PG_SubCortical_z_healthy_cont(:, 1:vertexnum, j) = (FLAIR_PG_SubCortical_temp(:, 1:vertexnum, case_cont_real(j)) - FLAIR_PG_SubCortical_mean) ./ FLAIR_PG_SubCortical_std;
                FLAIR_TG_SubCortical_z_healthy_cont(:, 1:vertexnum, j) = (FLAIR_TG_SubCortical_temp(:, 1:vertexnum, case_cont_real(j)) - FLAIR_TG_SubCortical_mean) ./ FLAIR_TG_SubCortical_std;
                FLAIR_RI_SubCortical_z_healthy_cont(isnan(FLAIR_RI_SubCortical_z_healthy_cont)) = 0; FLAIR_RI_SubCortical_z_healthy_cont(isinf(FLAIR_RI_SubCortical_z_healthy_cont)) = 0;
                FLAIR_PG_SubCortical_z_healthy_cont(isnan(FLAIR_PG_SubCortical_z_healthy_cont)) = 0; FLAIR_PG_SubCortical_z_healthy_cont(isinf(FLAIR_PG_SubCortical_z_healthy_cont)) = 0;
                FLAIR_TG_SubCortical_z_healthy_cont(isnan(FLAIR_TG_SubCortical_z_healthy_cont)) = 0; FLAIR_TG_SubCortical_z_healthy_cont(isinf(FLAIR_TG_SubCortical_z_healthy_cont)) = 0;
                fprintf([ 'RI, PG, TG in Sub\n' ]);
            end
            
        end
        
    end
    
    % disease controls
    for disease_control_computation = 1
        
        % read T1 files
        for T1_read_file = 1
            
            %% Intensity-based features
            T1_RI_IntCortical_temp2      = zeros(NumIntSurf+2, vertexnum, size(case_num_disease_cont, 1));
            T1_RI_SubCortical_temp2      = zeros(NumSubSurf, vertexnum, size(case_num_disease_cont, 1));
            T1_PG_IntCortical_temp2      = zeros(NumIntSurf+2, vertexnum, size(case_num_disease_cont, 1));
            T1_PG_gw_IntCortical_temp2   = zeros(size(case_num_disease_cont, 1), vertexnum);
            T1_PG_SubCortical_temp2      = zeros(NumSubSurf, vertexnum, size(case_num_disease_cont, 1));
            T1_TG_IntCortical_temp2      = zeros(NumIntSurf+2, vertexnum, size(case_num_disease_cont, 1));
            T1_TG_SubCortical_temp2      = zeros(NumSubSurf, vertexnum, size(case_num_disease_cont, 1));
            
            %% Morphological features
            T1_CT_midCortical_temp2      = zeros(size(case_num_disease_cont, 1), vertexnum);
            T1_SD_wmCortical_temp2       = zeros(size(case_num_disease_cont, 1), vertexnum);
            T1_MC_midCortical_temp2      = zeros(size(case_num_disease_cont, 1), vertexnum);
            
            data_dir = '/local_raid/seokjun/01_project/04_IntracorticalAnalysis/03_result/97_3T_TLE_HS/';
            file_postfix = [ '_quadratic_sm_5_rsl.txt' ];
            for j = 1 : size(case_num_disease_cont, 1)
                fprintf([ case_num_disease_cont{j} ' : ' ]);
                prefix_path = [ data_dir '/' case_num_disease_cont{j} '/measurement/' Prefix_disease_cont '_' case_num_disease_cont{j} ];
                postfix_surf = '_native_t1_';
                
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
                    
                    T1_RI_IntCortical_temp2(i, 1:vertexnum, j) = SurfStatReadData( { [ basename{1} 'RI_corrected' file_postfix ], [ basename{2} 'RI_corrected' file_postfix ] } );
                    T1_PG_IntCortical_temp2(i, 1:vertexnum, j) = SurfStatReadData( { [ basename{1} 'pg' file_postfix ], [ basename{2} 'pg' file_postfix ] } );
                    T1_TG_IntCortical_temp2(i, 1:vertexnum, j) = SurfStatReadData( { [ basename{1} 'tg' file_postfix ], [ basename{2} 'tg' file_postfix ] } );
                end
                fprintf([ 'T1: RI_corrected, PG, TG in Inc, ' ]);
                
                %% Morphological / Intensity-based features: MC, SD, PG_GW
                basename{1}  = [ prefix_path '_white_surface_left_'  num2str(NumMesh) postfix_surf ];
                basename{2}  = [ prefix_path '_white_surface_right_'  num2str(NumMesh) postfix_surf ];
                T1_PG_gw_IntCortical_temp2(j, 1:vertexnum) = SurfStatReadData( { [ basename{1} 'pg_GM_WM' file_postfix ], [ basename{2} 'pg_GM_WM' file_postfix ] } );
                
                basename{1}  = [ prefix_path '_intracortical_surface_2_left_'  num2str(NumMesh) postfix_surf ];
                basename{2}  = [ prefix_path '_intracortical_surface_2_right_'  num2str(NumMesh) postfix_surf ];
                T1_CT_midCortical_temp2(j, 1:vertexnum) = SurfStatReadData( { [ basename{1} 'ct' file_postfix ], [ basename{2} 'ct' file_postfix ] } );
                
                basename{1}  = [ prefix_path '_intracortical_surface_2_left_'  num2str(NumMesh) postfix_surf ];
                basename{2}  = [ prefix_path '_intracortical_surface_2_right_'  num2str(NumMesh) postfix_surf ];
                T1_MC_midCortical_temp2(j, 1:vertexnum) = SurfStatReadData( { [ basename{1} 'mc' file_postfix ], [ basename{2} 'mc' file_postfix ] } );
                
                basename{1}  = [ prefix_path '_white_surface_left_'  num2str(NumMesh) postfix_surf ];
                basename{2}  = [ prefix_path '_white_surface_right_'  num2str(NumMesh) postfix_surf ];
                T1_SD_wmCortical_temp2(j, 1:vertexnum) = SurfStatReadData( { [ basename{1} 'sd2' file_postfix ], [ basename{2} 'sd2' file_postfix ] } );
                fprintf([ 'CT, MC, SD, PG_gw in Inc, ' ]);
                
                %% Intensity-based features of subcortical surfaces: RI, pg, tg
                for i = 1 : NumSubSurf
                    basename{1}  = [ prefix_path '_white_surface_' num2str(i) '_left_'  num2str(NumMesh) postfix_surf ];
                    basename{2}  = [ prefix_path '_white_surface_' num2str(i) '_right_'  num2str(NumMesh) postfix_surf ];
                    
                    T1_RI_SubCortical_temp2(i, 1:vertexnum, j) = SurfStatReadData( { [ basename{1} 'RI' file_postfix ], [ basename{2} 'RI' file_postfix ] } );
                    T1_PG_SubCortical_temp2(i, 1:vertexnum, j) = SurfStatReadData( { [ basename{1} 'pg' file_postfix ], [ basename{2} 'pg' file_postfix ] } );
                    T1_TG_SubCortical_temp2(i, 1:vertexnum, j) = SurfStatReadData( { [ basename{1} 'tg' file_postfix ], [ basename{2} 'tg' file_postfix ] } );
                end
                fprintf([ 'RI, PG, TG in Sub\n' ]);
            end
            
        end

        % read FLAIR files
        for FLAIR_read_file = 1
            
            %% Intensity-based features
            FLAIR_RI_IntCortical_temp2      = zeros(NumIntSurf+2, vertexnum, size(case_num_disease_cont, 1));
            FLAIR_RI_SubCortical_temp2      = zeros(NumSubSurf, vertexnum, size(case_num_disease_cont, 1));
            FLAIR_PG_IntCortical_temp2      = zeros(NumIntSurf+2, vertexnum, size(case_num_disease_cont, 1));
            FLAIR_PG_gw_IntCortical_temp2   = zeros(size(case_num_disease_cont, 1), vertexnum);
            FLAIR_PG_SubCortical_temp2      = zeros(NumSubSurf, vertexnum, size(case_num_disease_cont, 1));
            FLAIR_TG_IntCortical_temp2      = zeros(NumIntSurf+2, vertexnum, size(case_num_disease_cont, 1));
            FLAIR_TG_SubCortical_temp2      = zeros(NumSubSurf, vertexnum, size(case_num_disease_cont, 1));
            
            for j = 1 : size(case_num_disease_cont, 1)
                fprintf([ case_num_disease_cont{j} ' : ' ]);
                prefix_path = [ data_dir '/' case_num_disease_cont{j} '/measurement/' Prefix_disease_cont '_' case_num_disease_cont{j} ];
                postfix_surf = '_native_flair_';
                
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
                    
                    FLAIR_RI_IntCortical_temp2(i, 1:vertexnum, j) = SurfStatReadData( { [ basename{1} 'RI_corrected' file_postfix ], [ basename{2} 'RI_corrected' file_postfix ] } );
                    FLAIR_PG_IntCortical_temp2(i, 1:vertexnum, j) = SurfStatReadData( { [ basename{1} 'pg' file_postfix ], [ basename{2} 'pg' file_postfix ] } );
                    FLAIR_TG_IntCortical_temp2(i, 1:vertexnum, j) = SurfStatReadData( { [ basename{1} 'tg' file_postfix ], [ basename{2} 'tg' file_postfix ] } );
                end
                fprintf([ 'FLAIR: RI_corrected, PG, TG in Inc, ' ]);
                
                %% Morphological / Intensity-based features: MC, SD, PG_GW
                basename{1}  = [ prefix_path '_white_surface_left_'  num2str(NumMesh) postfix_surf ];
                basename{2}  = [ prefix_path '_white_surface_right_'  num2str(NumMesh) postfix_surf ];
                FLAIR_PG_gw_IntCortical_temp2(j, 1:vertexnum) = SurfStatReadData( { [ basename{1} 'pg_GM_WM' file_postfix ], [ basename{2} 'pg_GM_WM' file_postfix ] } );
                
                fprintf([ 'PG_gw in Inc, ' ]);
                
                %% Intensity-based features of subcortical surfaces: RI, pg, tg
                for i = 1 : NumSubSurf
                    basename{1}  = [ prefix_path '_white_surface_' num2str(i) '_left_'  num2str(NumMesh) postfix_surf ];
                    basename{2}  = [ prefix_path '_white_surface_' num2str(i) '_right_'  num2str(NumMesh) postfix_surf ];
                    
                    FLAIR_RI_SubCortical_temp2(i, 1:vertexnum, j) = SurfStatReadData( { [ basename{1} 'RI' file_postfix ], [ basename{2} 'RI' file_postfix ] } );
                    FLAIR_PG_SubCortical_temp2(i, 1:vertexnum, j) = SurfStatReadData( { [ basename{1} 'pg' file_postfix ], [ basename{2} 'pg' file_postfix ] } );
                    FLAIR_TG_SubCortical_temp2(i, 1:vertexnum, j) = SurfStatReadData( { [ basename{1} 'tg' file_postfix ], [ basename{2} 'tg' file_postfix ] } );
                end
                fprintf([ 'RI, PG, TG in Sub\n' ]);
            end
            
        end

        % compute T1 z-score
        for T1_zscore_pat = 1
            
            %% z-score
            T1_RI_IntCortical_z_disease_cont         = zeros(NumIntSurf+2, vertexnum, size(case_num_disease_cont, 1));
            T1_RI_SubCortical_z_disease_cont         = zeros(NumSubSurf,   vertexnum, size(case_num_disease_cont, 1));
            T1_PG_IntCortical_z_disease_cont         = zeros(NumIntSurf+2, vertexnum, size(case_num_disease_cont, 1));
            T1_PG_SubCortical_z_disease_cont         = zeros(NumSubSurf,   vertexnum, size(case_num_disease_cont, 1));
            T1_TG_IntCortical_z_disease_cont         = zeros(NumIntSurf+2, vertexnum, size(case_num_disease_cont, 1));
            T1_TG_SubCortical_z_disease_cont         = zeros(NumSubSurf,   vertexnum, size(case_num_disease_cont, 1));
            T1_CT_midCortical_z_disease_cont         = zeros(size(case_num_disease_cont, 1), vertexnum);
            T1_MC_midCortical_z_disease_cont         = zeros(size(case_num_disease_cont, 1), vertexnum);
            T1_SD_wmCortical_z_disease_cont          = zeros(size(case_num_disease_cont, 1), vertexnum);
            T1_PG_gw_IntCortical_z_disease_cont      = zeros(size(case_num_disease_cont, 1), vertexnum);
            
            %% T1: calculate z-score of patients w.r.t control's distribution
            for j = 1 : size(case_num_disease_cont, 1)
                fprintf([ case_num_disease_cont{j} ' : ' ]);
                
                %% Intensity-based features: corrected RI, pg, tg
                T1_RI_IntCortical_z_disease_cont(:, 1:vertexnum, j) = (T1_RI_IntCortical_temp2(:, 1:vertexnum, j) - T1_RI_IntCortical_mean) ./ T1_RI_IntCortical_std;
                T1_PG_IntCortical_z_disease_cont(:, 1:vertexnum, j) = (T1_PG_IntCortical_temp2(:, 1:vertexnum, j) - T1_PG_IntCortical_mean) ./ T1_PG_IntCortical_std;
                T1_TG_IntCortical_z_disease_cont(:, 1:vertexnum, j) = (T1_TG_IntCortical_temp2(:, 1:vertexnum, j) - T1_TG_IntCortical_mean) ./ T1_TG_IntCortical_std;
                T1_RI_IntCortical_z_disease_cont(isnan(T1_RI_IntCortical_z_disease_cont)) = 0; T1_RI_IntCortical_z_disease_cont(isinf(T1_RI_IntCortical_z_disease_cont)) = 0;
                T1_PG_IntCortical_z_disease_cont(isnan(T1_PG_IntCortical_z_disease_cont)) = 0; T1_PG_IntCortical_z_disease_cont(isinf(T1_PG_IntCortical_z_disease_cont)) = 0;
                T1_TG_IntCortical_z_disease_cont(isnan(T1_TG_IntCortical_z_disease_cont)) = 0; T1_PG_IntCortical_z_disease_cont(isinf(T1_PG_IntCortical_z_disease_cont)) = 0;
                fprintf([ 'z-score of T1: RI_corrected, PG, TG in Inc, ' ]);
                
                %% Morphological / Intensity-based features: MC, SD, PG_GW
                T1_PG_gw_IntCortical_z_disease_cont(j, 1:vertexnum) = (T1_PG_gw_IntCortical_temp2(j, 1:vertexnum) - T1_PG_gw_IntCortical_mean) ./ T1_PG_gw_IntCortical_std;
                T1_CT_midCortical_z_disease_cont(j, 1:vertexnum)    = (T1_CT_midCortical_temp2(j, 1:vertexnum) - T1_CT_midCortical_mean) ./ T1_CT_midCortical_std;
                T1_MC_midCortical_z_disease_cont(j, 1:vertexnum)    = (T1_MC_midCortical_temp2(j, 1:vertexnum) - T1_MC_midCortical_mean) ./ T1_MC_midCortical_std;
                T1_SD_wmCortical_z_disease_cont(j, 1:vertexnum)     = (T1_SD_wmCortical_temp2(j, 1:vertexnum) - T1_SD_wmCortical_mean) ./ T1_SD_wmCortical_std;
                T1_PG_gw_IntCortical_z_disease_cont(isnan(T1_PG_gw_IntCortical_z_disease_cont)) = 0; T1_PG_gw_IntCortical_z_disease_cont(isinf(T1_PG_gw_IntCortical_z_disease_cont)) = 0;
                T1_CT_midCortical_z_disease_cont(isnan(T1_CT_midCortical_z_disease_cont)) = 0; T1_CT_midCortical_z_disease_cont(isinf(T1_CT_midCortical_z_disease_cont)) = 0;
                T1_MC_midCortical_z_disease_cont(isnan(T1_MC_midCortical_z_disease_cont)) = 0; T1_MC_midCortical_z_disease_cont(isinf(T1_MC_midCortical_z_disease_cont)) = 0;
                T1_SD_wmCortical_z_disease_cont(isnan(T1_SD_wmCortical_z_disease_cont)) = 0; T1_SD_wmCortical_z_disease_cont(isinf(T1_SD_wmCortical_z_disease_cont)) = 0;
                fprintf([ 'CT, MC, SD, PG_gw in Inc, ' ]);
                
                %% Intensity-based features of subcortical surfaces: RI, pg, tg
                T1_RI_SubCortical_z_disease_cont(:, 1:vertexnum, j) = (T1_RI_SubCortical_temp2(:, 1:vertexnum, j) - T1_RI_SubCortical_mean) ./ T1_RI_SubCortical_std;
                T1_PG_SubCortical_z_disease_cont(:, 1:vertexnum, j) = (T1_PG_SubCortical_temp2(:, 1:vertexnum, j) - T1_PG_SubCortical_mean) ./ T1_PG_SubCortical_std;
                T1_TG_SubCortical_z_disease_cont(:, 1:vertexnum, j) = (T1_TG_SubCortical_temp2(:, 1:vertexnum, j) - T1_TG_SubCortical_mean) ./ T1_TG_SubCortical_std;
                T1_RI_SubCortical_z_disease_cont(isnan(T1_RI_SubCortical_z_disease_cont)) = 0; T1_RI_SubCortical_z_disease_cont(isinf(T1_RI_SubCortical_z_disease_cont)) = 0;
                T1_PG_SubCortical_z_disease_cont(isnan(T1_PG_SubCortical_z_disease_cont)) = 0; T1_PG_SubCortical_z_disease_cont(isinf(T1_PG_SubCortical_z_disease_cont)) = 0;
                T1_TG_SubCortical_z_disease_cont(isnan(T1_TG_SubCortical_z_disease_cont)) = 0; T1_TG_SubCortical_z_disease_cont(isinf(T1_TG_SubCortical_z_disease_cont)) = 0;
                fprintf([ 'RI, PG, TG in Sub\n' ]);
            end
            
        end

        % compute FLAIR z-score
        for FLAIR_z_disease_contscore_pat = 1
            
            FLAIR_RI_IntCortical_z_disease_cont         = zeros(NumIntSurf+2, vertexnum, size(case_num_disease_cont, 1));
            FLAIR_RI_SubCortical_z_disease_cont         = zeros(NumSubSurf,   vertexnum, size(case_num_disease_cont, 1));
            FLAIR_PG_IntCortical_z_disease_cont         = zeros(NumIntSurf+2, vertexnum, size(case_num_disease_cont, 1));
            FLAIR_PG_SubCortical_z_disease_cont         = zeros(NumSubSurf,   vertexnum, size(case_num_disease_cont, 1));
            FLAIR_TG_IntCortical_z_disease_cont         = zeros(NumIntSurf+2, vertexnum, size(case_num_disease_cont, 1));
            FLAIR_TG_SubCortical_z_disease_cont         = zeros(NumSubSurf,   vertexnum, size(case_num_disease_cont, 1));
            FLAIR_PG_gw_IntCortical_z_disease_cont      = zeros(size(case_num_disease_cont, 1), vertexnum);
            %% FLAIR: calculate z-score of patients w.r.t control's distribution
            for j = 1 : size(case_num_disease_cont, 1)
                fprintf([ case_num_disease_cont{j} ' : ' ]);
                
                %% Intensity-based features: corrected RI, pg, tg
                FLAIR_RI_IntCortical_z_disease_cont(:, 1:vertexnum, j) = (FLAIR_RI_IntCortical_temp2(:, 1:vertexnum, j) - FLAIR_RI_IntCortical_mean) ./ FLAIR_RI_IntCortical_std;
                FLAIR_PG_IntCortical_z_disease_cont(:, 1:vertexnum, j) = (FLAIR_PG_IntCortical_temp2(:, 1:vertexnum, j) - FLAIR_PG_IntCortical_mean) ./ FLAIR_PG_IntCortical_std;
                FLAIR_TG_IntCortical_z_disease_cont(:, 1:vertexnum, j) = (FLAIR_TG_IntCortical_temp2(:, 1:vertexnum, j) - FLAIR_TG_IntCortical_mean) ./ FLAIR_TG_IntCortical_std;
                FLAIR_RI_IntCortical_z_disease_cont(isnan(FLAIR_RI_IntCortical_z_disease_cont)) = 0; FLAIR_RI_IntCortical_z_disease_cont(isinf(FLAIR_RI_IntCortical_z_disease_cont)) = 0;
                FLAIR_PG_IntCortical_z_disease_cont(isnan(FLAIR_PG_IntCortical_z_disease_cont)) = 0; FLAIR_PG_IntCortical_z_disease_cont(isinf(FLAIR_PG_IntCortical_z_disease_cont)) = 0;
                FLAIR_TG_IntCortical_z_disease_cont(isnan(FLAIR_TG_IntCortical_z_disease_cont)) = 0; FLAIR_TG_IntCortical_z_disease_cont(isinf(FLAIR_TG_IntCortical_z_disease_cont)) = 0;
                fprintf([ 'z-score of FLAIR: RI_corrected, PG, TG in Inc, ' ]);
                
                %% Morphological / Intensity-based features: MC, SD, PG_GW
                FLAIR_PG_gw_IntCortical_z_disease_cont(j, 1:vertexnum) = (FLAIR_PG_gw_IntCortical_temp2(j, 1:vertexnum) - FLAIR_PG_gw_IntCortical_mean) ./ FLAIR_PG_gw_IntCortical_std;
                FLAIR_PG_gw_IntCortical_z_disease_cont(isnan(FLAIR_PG_gw_IntCortical_z_disease_cont)) = 0;
                fprintf([ 'PG_gw in Inc, ' ]);
                
                %% Intensity-based features of subcortical surfaces: RI, pg, tg
                FLAIR_RI_SubCortical_z_disease_cont(:, 1:vertexnum, j) = (FLAIR_RI_SubCortical_temp2(:, 1:vertexnum, j) - FLAIR_RI_SubCortical_mean) ./ FLAIR_RI_SubCortical_std;
                FLAIR_PG_SubCortical_z_disease_cont(:, 1:vertexnum, j) = (FLAIR_PG_SubCortical_temp2(:, 1:vertexnum, j) - FLAIR_PG_SubCortical_mean) ./ FLAIR_PG_SubCortical_std;
                FLAIR_TG_SubCortical_z_disease_cont(:, 1:vertexnum, j) = (FLAIR_TG_SubCortical_temp2(:, 1:vertexnum, j) - FLAIR_TG_SubCortical_mean) ./ FLAIR_TG_SubCortical_std;
                FLAIR_RI_SubCortical_z_disease_cont(isnan(FLAIR_RI_SubCortical_z_disease_cont)) = 0; FLAIR_RI_SubCortical_z_disease_cont(isinf(FLAIR_RI_SubCortical_z_disease_cont)) = 0;
                FLAIR_PG_SubCortical_z_disease_cont(isnan(FLAIR_PG_SubCortical_z_disease_cont)) = 0; FLAIR_PG_SubCortical_z_disease_cont(isinf(FLAIR_PG_SubCortical_z_disease_cont)) = 0;
                FLAIR_TG_SubCortical_z_disease_cont(isnan(FLAIR_TG_SubCortical_z_disease_cont)) = 0; FLAIR_TG_SubCortical_z_disease_cont(isinf(FLAIR_TG_SubCortical_z_disease_cont)) = 0;
                fprintf([ 'RI, PG, TG in Sub\n' ]);
            end
            
        end
        
        for visualization_zscore = 1
            
            %% 1st mean and SD examination ...
            IntMask = repmat(TemplateMask, [5, 1]);
            SubMask = repmat(TemplateMask, [3, 1]);
            OUTPATH_temp = '/local_raid/seokjun/01_project/04_IntracorticalAnalysis/03_result/zscoremap/';
            for j = 1 : size(case_num_disease_cont, 1)
                
                %     f = figure; CSFSurfStatView(T1_RI_IntCortical_z_disease_cont(:, 1:vertexnum, idx).*IntMask, TemplateSurf(1:5), 'INC', 'RI z-score');         CSFSurfStatViewColLim([-5 5]);
                %     exportfigbo(f,[OUTPATH_temp '/' case_num_disease_cont{j} '/' Prefix_disease_cont '_' case_num_disease_cont{j} '_z_disease_cont-score_T1_Intra_RI.png' ], 'png', 6); close(f);
                %     f = figure; CSFSurfStatView(T1_PG_IntCortical_z_disease_cont(:, 1:vertexnum, idx).*IntMask, TemplateSurf(1:5), 'INC', 'PG z-score');         CSFSurfStatViewColLim([-5 5]);
                %     exportfigbo(f,[OUTPATH_temp '/' case_num_disease_cont{j} '/' Prefix_disease_cont '_' case_num_disease_cont{j} '_z_disease_cont-score_T1_Intra_PG.png' ], 'png', 6); close(f);
                %     f = figure; CSFSurfStatView(T1_TG_IntCortical_z_disease_cont(:, 1:vertexnum, idx).*IntMask, TemplateSurf(1:5), 'INC', 'TG z-score');         CSFSurfStatViewColLim([-5 5]);
                %     exportfigbo(f,[OUTPATH_temp '/' case_num_disease_cont{j} '/' Prefix_disease_cont '_' case_num_disease_cont{j} '_z_disease_cont-score_T1_Intra_TG.png' ], 'png', 6); close(f);
                %     f = figure; CSFSurfStatView(T1_PG_gw_IntCortical_z_disease_cont(idx, 1:vertexnum).*TemplateMask, TemplateSurf{3}, 'MF', 'PG z-score');     CSFSurfStatViewColLim([-5 5]);
                %     exportfigbo(f,[OUTPATH_temp '/' case_num_disease_cont{j} '/' Prefix_disease_cont '_' case_num_disease_cont{j} '_z_disease_cont-score_T1_GM-WM_PG.png' ], 'png', 6); close(f);
                %     f = figure; CSFSurfStatView(T1_CT_midCortical_z_disease_cont(idx, 1:vertexnum).*TemplateMask, TemplateSurf{3}, 'MF', 'CT z-score');          CSFSurfStatViewColLim([-5 5]);
                %     exportfigbo(f,[OUTPATH_temp '/' case_num_disease_cont{j} '/' Prefix_disease_cont '_' case_num_disease_cont{j} '_z_disease_cont-score_T1_Mid_CT.png' ], 'png', 6); close(f);
                %     f = figure; CSFSurfStatView(T1_MC_midCortical_z_disease_cont(idx, 1:vertexnum).*TemplateMask, TemplateSurf{3}, 'MF', 'MC z-score');          CSFSurfStatViewColLim([-5 5]);
                %     exportfigbo(f,[OUTPATH_temp '/' case_num_disease_cont{j} '/' Prefix_disease_cont '_' case_num_disease_cont{j} '_z_disease_cont-score_T1_Mid_MC.png' ], 'png', 6); close(f);
                %     f = figure; CSFSurfStatView(T1_SD_wmCortical_z_disease_cont(idx, 1:vertexnum).*TemplateMask, TemplateSurf{3},  'MF', 'SD z-score');          CSFSurfStatViewColLim([-5 5]);
                %     exportfigbo(f,[OUTPATH_temp '/' case_num_disease_cont{j} '/' Prefix_disease_cont '_' case_num_disease_cont{j} '_z_disease_cont-score_T1_WM_SD.png' ], 'png', 6); close(f);
                %
                %     f = figure; CSFSurfStatView(T1_RI_SubCortical_z_disease_cont(:, 1:vertexnum, idx).*SubMask, TemplateSurf(6:8), 'SBC', 'RI z-score');         CSFSurfStatViewColLim([-5 5]);
                %     exportfigbo(f,[OUTPATH_temp '/' case_num_disease_cont{j} '/' Prefix_disease_cont '_' case_num_disease_cont{j} '_z_disease_cont-score_T1_Sub_RI.png' ], 'png', 6); close(f);
                %     f = figure; CSFSurfStatView(T1_PG_SubCortical_z_disease_cont(:, 1:vertexnum, idx).*SubMask, TemplateSurf(6:8), 'SBC', 'PG z-score');         CSFSurfStatViewColLim([-5 5]);
                %     exportfigbo(f,[OUTPATH_temp '/' case_num_disease_cont{j} '/' Prefix_disease_cont '_' case_num_disease_cont{j} '_z_disease_cont-score_T1_Sub_PG.png' ], 'png', 6); close(f);
                %     f = figure; CSFSurfStatView(T1_TG_SubCortical_z_disease_cont(:, 1:vertexnum, idx).*SubMask, TemplateSurf(6:8), 'SBC', 'TG z-score');         CSFSurfStatViewColLim([-5 5]);
                %     exportfigbo(f,[OUTPATH_temp '/' case_num_disease_cont{j} '/' Prefix_disease_cont '_' case_num_disease_cont{j} '_z_disease_cont-score_T1_Sub_TG.png' ], 'png', 6); close(f);
                %
                %     f = figure; CSFSurfStatView(FLAIR_RI_IntCortical_z_disease_cont(:, 1:vertexnum, idx).*IntMask, TemplateSurf(1:5), 'INC', 'RI z-score');      CSFSurfStatViewColLim([-5 5]);
                %     exportfigbo(f,[OUTPATH_temp '/' case_num_disease_cont{j} '/' Prefix_disease_cont '_' case_num_disease_cont{j} '_z_disease_cont-score_FLAIR_Intra_RI.png' ], 'png', 6); close(f);
                %     f = figure; CSFSurfStatView(FLAIR_PG_IntCortical_z_disease_cont(:, 1:vertexnum, idx).*IntMask, TemplateSurf(1:5), 'INC', 'PG z-score');      CSFSurfStatViewColLim([-5 5]);
                %     exportfigbo(f,[OUTPATH_temp '/' case_num_disease_cont{j} '/' Prefix_disease_cont '_' case_num_disease_cont{j} '_z_disease_cont-score_FLAIR_Intra_PG.png' ], 'png', 6); close(f);
                %     f = figure; CSFSurfStatView(FLAIR_TG_IntCortical_z_disease_cont(:, 1:vertexnum, idx).*IntMask, TemplateSurf(1:5), 'INC', 'TG z-score');      CSFSurfStatViewColLim([-5 5]);
                %     exportfigbo(f,[OUTPATH_temp '/' case_num_disease_cont{j} '/' Prefix_disease_cont '_' case_num_disease_cont{j} '_z_disease_cont-score_FLAIR_Intra_TG.png' ], 'png', 6); close(f);
                %     f = figure; CSFSurfStatView(FLAIR_PG_gw_IntCortical_z_disease_cont(idx, 1:vertexnum).*TemplateMask, TemplateSurf{3}, 'MF', 'PG z-score'); CSFSurfStatViewColLim([-5 5]);
                %     exportfigbo(f,[OUTPATH_temp '/' case_num_disease_cont{j} '/' Prefix_disease_cont '_' case_num_disease_cont{j} '_z_disease_cont-score_FLAIR_GM-WM_PG.png' ], 'png', 6); close(f);
                %
                %     f = figure; CSFSurfStatView(FLAIR_RI_SubCortical_z_disease_cont(:, 1:vertexnum, idx).*SubMask, TemplateSurf(6:8), 'SBC', 'RI z-score');      CSFSurfStatViewColLim([-5 5]);
                %     exportfigbo(f,[OUTPATH_temp '/' case_num_disease_cont{j} '/' Prefix_disease_cont '_' case_num_disease_cont{j} '_z_disease_cont-score_FLAIR_Sub_RI.png' ], 'png', 6); close(f);
                %     f = figure; CSFSurfStatView(FLAIR_PG_SubCortical_z_disease_cont(:, 1:vertexnum, idx).*SubMask, TemplateSurf(6:8), 'SBC', 'PG z-score');      CSFSurfStatViewColLim([-5 5]);
                %     exportfigbo(f,[OUTPATH_temp '/' case_num_disease_cont{j} '/' Prefix_disease_cont '_' case_num_disease_cont{j} '_z_disease_cont-score_FLAIR_Sub_PG.png' ], 'png', 6); close(f);
                %     f = figure; CSFSurfStatView(FLAIR_TG_SubCortical_z_disease_cont(:, 1:vertexnum, idx).*SubMask, TemplateSurf(6:8), 'SBC', 'TG z-score');      CSFSurfStatViewColLim([-5 5]);
                %     exportfigbo(f,[OUTPATH_temp '/' case_num_disease_cont{j} '/' Prefix_disease_cont '_' case_num_disease_cont{j} '_z_disease_cont-score_FLAIR_Sub_TG.png' ], 'png', 6); close(f);
                
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
        mean_z_lesion_disease_cont = zeros(53, size(case_num_disease_cont, 1));
        std_z_lesion_disease_cont  = zeros(53, size(case_num_disease_cont, 1));
        case_num_pat_temp = case_num_pat;
        for j = 1 : size(case_num_pat_temp, 1)
            idx = find(strcmp(case_num_pat_temp, case_num_pat_temp{j}));
            
            if(exist([ lesion_label_dir '/' Prefix_pat '_' case_num_pat_temp{j} '_label_union_left_rsl.txt' ], 'file'))
                lesion_label_data = SurfStatReadData( { [ lesion_label_dir '/' Prefix_pat '_' case_num_pat_temp{j} '_label_union_left_rsl.txt' ], ...
                                                        [ lesion_label_dir '/' Prefix_pat '_' case_num_pat_temp{j} '_label_union_right_rsl.txt' ] } );
                
                lesion_label_v = find(lesion_label_data ~= 0);
                
                for c = 1 : size(case_num_disease_cont, 1)
                    %% Cortical T1 intensity
                    mean_z_lesion_disease_cont(1:5, j, c) = mean(T1_RI_IntCortical_z_disease_cont(:, lesion_label_v, c), 2);
                    std_z_lesion_disease_cont(1:5, j, c)  = std(T1_RI_IntCortical_z_disease_cont(:, lesion_label_v, c), 0, 2);
                    mean_z_lesion_disease_cont(6:10, j, c) = mean(T1_PG_IntCortical_z_disease_cont(:, lesion_label_v, c), 2);
                    std_z_lesion_disease_cont(6:10, j, c)  = std(T1_PG_IntCortical_z_disease_cont(:, lesion_label_v, c), 0, 2);
                    mean_z_lesion_disease_cont(11:15, j, c) = mean(T1_TG_IntCortical_z_disease_cont(:, lesion_label_v, c), 2);
                    std_z_lesion_disease_cont(11:15, j, c)  = std(T1_TG_IntCortical_z_disease_cont(:, lesion_label_v, c), 0, 2);
                    mean_z_lesion_disease_cont(16, j, c) = mean(T1_PG_gw_IntCortical_z_disease_cont(c, lesion_label_v), 2);
                    std_z_lesion_disease_cont(16, j, c)  = std(T1_PG_gw_IntCortical_z_disease_cont(c, lesion_label_v), 0, 2);
                    
                    %% T1 Morphological metrics
                    mean_z_lesion_disease_cont(17, j, c) = mean(T1_CT_midCortical_z_disease_cont(c, lesion_label_v), 2);
                    std_z_lesion_disease_cont(17, j, c)  = std(T1_CT_midCortical_z_disease_cont(c, lesion_label_v), 0, 2);
                    mean_z_lesion_disease_cont(18, j, c) = mean(T1_MC_midCortical_z_disease_cont(c, lesion_label_v), 2);
                    std_z_lesion_disease_cont(18, j, c)  = std(T1_MC_midCortical_z_disease_cont(c, lesion_label_v), 0, 2);
                    mean_z_lesion_disease_cont(19, j, c) = mean(T1_SD_wmCortical_z_disease_cont(c, lesion_label_v), 2);
                    std_z_lesion_disease_cont(19, j, c)  = std(T1_SD_wmCortical_z_disease_cont(c, lesion_label_v), 0, 2);
                    
                    %% subcortical T1 intensity
                    mean_z_lesion_disease_cont(20:22, j, c) = mean(T1_RI_SubCortical_z_disease_cont(:, lesion_label_v, c), 2);
                    std_z_lesion_disease_cont(20:22, j, c)  = std(T1_RI_SubCortical_z_disease_cont(:, lesion_label_v, c), 0, 2);
                    mean_z_lesion_disease_cont(23:25, j, c) = mean(T1_PG_SubCortical_z_disease_cont(:, lesion_label_v, c), 2);
                    std_z_lesion_disease_cont(23:25, j, c)  = std(T1_PG_SubCortical_z_disease_cont(:, lesion_label_v, c), 0, 2);
                    mean_z_lesion_disease_cont(26:28, j, c) = mean(T1_TG_SubCortical_z_disease_cont(:, lesion_label_v, c), 2);
                    std_z_lesion_disease_cont(26:28, j, c)  = std(T1_TG_SubCortical_z_disease_cont(:, lesion_label_v, c), 0, 2);
                    
                    %% Cortical FLAIR intensity
                    mean_z_lesion_disease_cont(29:33, j, c) = mean(FLAIR_RI_IntCortical_z_disease_cont(:, lesion_label_v, c), 2);
                    std_z_lesion_disease_cont(29:33, j, c)  = std(FLAIR_RI_IntCortical_z_disease_cont(:, lesion_label_v, c), 0, 2);
                    mean_z_lesion_disease_cont(34:38, j, c) = mean(FLAIR_PG_IntCortical_z_disease_cont(:, lesion_label_v, c), 2);
                    std_z_lesion_disease_cont(34:38, j, c)  = std(FLAIR_PG_IntCortical_z_disease_cont(:, lesion_label_v, c), 0, 2);
                    mean_z_lesion_disease_cont(39:43, j, c) = mean(FLAIR_TG_IntCortical_z_disease_cont(:, lesion_label_v, c), 2);
                    std_z_lesion_disease_cont(39:43, j, c)  = std(FLAIR_TG_IntCortical_z_disease_cont(:, lesion_label_v, c), 0, 2);
                    mean_z_lesion_disease_cont(44, j, c) = mean(FLAIR_PG_gw_IntCortical_z_disease_cont(c, lesion_label_v), 2);
                    std_z_lesion_disease_cont(44, j, c)  = std(FLAIR_PG_gw_IntCortical_z_disease_cont(c, lesion_label_v), 0, 2);
                    
                    %% subcortical FLAIR intensity
                    mean_z_lesion_disease_cont(45:47, j, c) = mean(FLAIR_RI_SubCortical_z_disease_cont(:, lesion_label_v, c), 2);
                    std_z_lesion_disease_cont(45:47, j, c)  = std(FLAIR_RI_SubCortical_z_disease_cont(:, lesion_label_v, c), 0, 2);
                    mean_z_lesion_disease_cont(48:50, j, c) = mean(FLAIR_PG_SubCortical_z_disease_cont(:, lesion_label_v, c), 2);
                    std_z_lesion_disease_cont(48:50, j, c)  = std(FLAIR_PG_SubCortical_z_disease_cont(:, lesion_label_v, c), 0, 2);
                    mean_z_lesion_disease_cont(51:53, j, c) = mean(FLAIR_TG_SubCortical_z_disease_cont(:, lesion_label_v, c), 2);
                    std_z_lesion_disease_cont(51:53, j, c)  = std(FLAIR_TG_SubCortical_z_disease_cont(:, lesion_label_v, c), 0, 2);
                end
            else
                disp([ case_num_pat_temp{j} ' does not have a lesion label file' ]);
            end
        end
        
        % healthy controls        
        mean_z_lesion_healthy_cont = zeros(53, size(case_num_healthy_cont, 1), size(case_num_cont_temp, 1));
        std_z_lesion_healthy_cont  = zeros(53, size(case_num_healthy_cont, 1), size(case_num_cont_temp, 1));        
        for j = 1 : size(case_num_pat_temp, 1)
            idx = find(strcmp(case_num_pat_temp, case_num_pat_temp{j}));
            
            if(exist([ lesion_label_dir '/' Prefix_pat '_' case_num_pat_temp{j} '_label_union_left_rsl.txt' ], 'file'))
                lesion_label_data = SurfStatReadData( { [ lesion_label_dir '/' Prefix_pat '_' case_num_pat_temp{j} '_label_union_left_rsl.txt' ], ...
                                                        [ lesion_label_dir '/' Prefix_pat '_' case_num_pat_temp{j} '_label_union_right_rsl.txt' ] } );
                
                lesion_label_v = find(lesion_label_data ~= 0);
                
                for c = 1 : size(case_num_cont_temp, 1)
                    %% Cortical T1 intensity
                    mean_z_lesion_healthy_cont(1:5, j, c) = mean(T1_RI_IntCortical_z_healthy_cont(:, lesion_label_v, c), 2);
                    std_z_lesion_healthy_cont(1:5, j, c)  = std(T1_RI_IntCortical_z_healthy_cont(:, lesion_label_v, c), 0, 2);
                    mean_z_lesion_healthy_cont(6:10, j, c) = mean(T1_PG_IntCortical_z_healthy_cont(:, lesion_label_v, c), 2);
                    std_z_lesion_healthy_cont(6:10, j, c)  = std(T1_PG_IntCortical_z_healthy_cont(:, lesion_label_v, c), 0, 2);
                    mean_z_lesion_healthy_cont(11:15, j, c) = mean(T1_TG_IntCortical_z_healthy_cont(:, lesion_label_v, c), 2);
                    std_z_lesion_healthy_cont(11:15, j, c)  = std(T1_TG_IntCortical_z_healthy_cont(:, lesion_label_v, c), 0, 2);
                    mean_z_lesion_healthy_cont(16, j, c) = mean(T1_PG_gw_IntCortical_z_healthy_cont(c, lesion_label_v), 2);
                    std_z_lesion_healthy_cont(16, j, c)  = std(T1_PG_gw_IntCortical_z_healthy_cont(c, lesion_label_v), 0, 2);
                    
                    %% T1 Morphological metrics
                    mean_z_lesion_healthy_cont(17, j, c) = mean(T1_CT_midCortical_z_healthy_cont(c, lesion_label_v), 2);
                    std_z_lesion_healthy_cont(17, j, c)  = std(T1_CT_midCortical_z_healthy_cont(c, lesion_label_v), 0, 2);
                    mean_z_lesion_healthy_cont(18, j, c) = mean(T1_MC_midCortical_z_healthy_cont(c, lesion_label_v), 2);
                    std_z_lesion_healthy_cont(18, j, c)  = std(T1_MC_midCortical_z_healthy_cont(c, lesion_label_v), 0, 2);
                    mean_z_lesion_healthy_cont(19, j, c) = mean(T1_SD_wmCortical_z_healthy_cont(c, lesion_label_v), 2);
                    std_z_lesion_healthy_cont(19, j, c)  = std(T1_SD_wmCortical_z_healthy_cont(c, lesion_label_v), 0, 2);
                    
                    %% subcortical T1 intensity
                    mean_z_lesion_healthy_cont(20:22, j, c) = mean(T1_RI_SubCortical_z_healthy_cont(:, lesion_label_v, c), 2);
                    std_z_lesion_healthy_cont(20:22, j, c)  = std(T1_RI_SubCortical_z_healthy_cont(:, lesion_label_v, c), 0, 2);
                    mean_z_lesion_healthy_cont(23:25, j, c) = mean(T1_PG_SubCortical_z_healthy_cont(:, lesion_label_v, c), 2);
                    std_z_lesion_healthy_cont(23:25, j, c)  = std(T1_PG_SubCortical_z_healthy_cont(:, lesion_label_v, c), 0, 2);
                    mean_z_lesion_healthy_cont(26:28, j, c) = mean(T1_TG_SubCortical_z_healthy_cont(:, lesion_label_v, c), 2);
                    std_z_lesion_healthy_cont(26:28, j, c)  = std(T1_TG_SubCortical_z_healthy_cont(:, lesion_label_v, c), 0, 2);
                    
                    %% Cortical FLAIR intensity
                    mean_z_lesion_healthy_cont(29:33, j, c) = mean(FLAIR_RI_IntCortical_z_healthy_cont(:, lesion_label_v, c), 2);
                    std_z_lesion_healthy_cont(29:33, j, c)  = std(FLAIR_RI_IntCortical_z_healthy_cont(:, lesion_label_v, c), 0, 2);
                    mean_z_lesion_healthy_cont(34:38, j, c) = mean(FLAIR_PG_IntCortical_z_healthy_cont(:, lesion_label_v, c), 2);
                    std_z_lesion_healthy_cont(34:38, j, c)  = std(FLAIR_PG_IntCortical_z_healthy_cont(:, lesion_label_v, c), 0, 2);
                    mean_z_lesion_healthy_cont(39:43, j, c) = mean(FLAIR_TG_IntCortical_z_healthy_cont(:, lesion_label_v, c), 2);
                    std_z_lesion_healthy_cont(39:43, j, c)  = std(FLAIR_TG_IntCortical_z_healthy_cont(:, lesion_label_v, c), 0, 2);
                    mean_z_lesion_healthy_cont(44, j, c) = mean(FLAIR_PG_gw_IntCortical_z_healthy_cont(c, lesion_label_v), 2);
                    std_z_lesion_healthy_cont(44, j, c)  = std(FLAIR_PG_gw_IntCortical_z_healthy_cont(c, lesion_label_v), 0, 2);
                    
                    %% subcortical FLAIR intensity
                    mean_z_lesion_healthy_cont(45:47, j, c) = mean(FLAIR_RI_SubCortical_z_healthy_cont(:, lesion_label_v, c), 2);
                    std_z_lesion_healthy_cont(45:47, j, c)  = std(FLAIR_RI_SubCortical_z_healthy_cont(:, lesion_label_v, c), 0, 2);
                    mean_z_lesion_healthy_cont(48:50, j, c) = mean(FLAIR_PG_SubCortical_z_healthy_cont(:, lesion_label_v, c), 2);
                    std_z_lesion_healthy_cont(48:50, j, c)  = std(FLAIR_PG_SubCortical_z_healthy_cont(:, lesion_label_v, c), 0, 2);
                    mean_z_lesion_healthy_cont(51:53, j, c) = mean(FLAIR_TG_SubCortical_z_healthy_cont(:, lesion_label_v, c), 2);
                    std_z_lesion_healthy_cont(51:53, j, c)  = std(FLAIR_TG_SubCortical_z_healthy_cont(:, lesion_label_v, c), 0, 2);
                end
            else
                disp([ case_num_pat_temp{j} ' does not have a lesion label file' ]);
            end
        end
        
    end
    
    for compute_lesion_size = 1
        
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
        
        if(~exist('zscore_database_T1_FLAIR_disease_cont.mat', 'file'))
            save('zscore_database_T1_FLAIR_disease_cont.mat', ...
                'T1_RI_IntCortical_z_disease_cont', 'T1_RI_SubCortical_z_disease_cont', 'T1_PG_IntCortical_z_disease_cont', 'T1_PG_SubCortical_z_disease_cont', 'T1_TG_IntCortical_z_disease_cont', ...
                'T1_TG_SubCortical_z_disease_cont', 'T1_CT_midCortical_z_disease_cont', 'T1_MC_midCortical_z_disease_cont', 'T1_SD_wmCortical_z_disease_cont', 'T1_PG_gw_IntCortical_z_disease_cont', ...
                'FLAIR_RI_IntCortical_z_disease_cont', 'FLAIR_RI_SubCortical_z_disease_cont', 'FLAIR_PG_IntCortical_z_disease_cont', 'FLAIR_PG_SubCortical_z_disease_cont', ...
                'FLAIR_TG_IntCortical_z_disease_cont', 'FLAIR_TG_SubCortical_z_disease_cont', 'FLAIR_PG_gw_IntCortical_z_disease_cont', ...
                'mean_z_lesion_disease_cont', 'std_z_lesion_disease_cont');
                
        end
        
        if(~exist('zscore_database_T1_FLAIR_healthy_cont.mat', 'file'))
            save('zscore_database_T1_FLAIR_healthy_cont.mat', ...
                'T1_RI_IntCortical_z_healthy_cont', 'T1_RI_SubCortical_z_healthy_cont', 'T1_PG_IntCortical_z_healthy_cont', 'T1_PG_SubCortical_z_healthy_cont', 'T1_TG_IntCortical_z_healthy_cont', ...
                'T1_TG_SubCortical_z_healthy_cont', 'T1_CT_midCortical_z_healthy_cont', 'T1_MC_midCortical_z_healthy_cont', 'T1_SD_wmCortical_z_healthy_cont', 'T1_PG_gw_IntCortical_z_healthy_cont', ...
                'FLAIR_RI_IntCortical_z_healthy_cont', 'FLAIR_RI_SubCortical_z_healthy_cont', 'FLAIR_PG_IntCortical_z_healthy_cont', 'FLAIR_PG_SubCortical_z_healthy_cont', ...
                'FLAIR_TG_IntCortical_z_healthy_cont', 'FLAIR_TG_SubCortical_z_healthy_cont', 'FLAIR_PG_gw_IntCortical_z_healthy_cont', ...
                'mean_z_lesion_healthy_cont', 'std_z_lesion_healthy_cont');
        end
        
    end

end
