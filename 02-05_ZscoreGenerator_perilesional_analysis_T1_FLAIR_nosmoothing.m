clear all
close all
addpath('/data/noel/noel2/local/matlab/surfstat_chicago');

%% Read demograpic data
for read_demodata = 1
    
    OutDir                    = '/local_raid/seokjun/01_project/04_IntracorticalAnalysis/03_result/'
    Group_cont                = 'control'
    Prefix_cont               = 'TLE'
    Cases_pat                 = '/local_raid/seokjun/01_project/04_IntracorticalAnalysis/01_analysis/Demographic_data_FCD.txt'
    Group_pat                 = 'FCD'
    Prefix_pat                = 'mcd'
    Cases_cont                = '/local_raid/seokjun/01_project/04_IntracorticalAnalysis/01_analysis/Demographic_data_control.txt'
    
    Left_Right                = 'both'
    NumIntSurf                = 3
    NumSubSurf                = 3
    NumMesh                   = 81920
    SamplingSpace             = 'native'
    img_contrast              = {'t1'}
    Kernel                    = 1
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
    
    % Modality: T1, FLAIR, IR
    
    %% Feature: RI_IntraCortical (vertexnum x 5 x numControl),
    %%          RI_SubCortical   (vertexnum x 3 x numControl),
    %%          PG_IntraCortical (vertexnum x 5 x numControl),
    %%          PG_GM_WM         (vertexnum x numControl),
    %%          PG_SubCortical   (vertexnum x 3 x numControl),
    %%          TG_IntraCortical (vertexnum x 5 x numControl),
    %%          TG_SubCortical   (vertexnum x 3 x numControl),
    %%          CT, MC, SD       (vertexnum x numControl, only T1)
    
    %% Abb: RI relative intensity
    %%      PG perpendicular gradient
    %%      TG tangential gradient
    %%      CT cortical thickness
    %%      MC mean curvature
    %%      SD sulcal depth

    % controls
    for control_computation = 1
        
        % T1
        for T1_read_file = 1
            
            T1_RI_IntCortical_temp      = zeros(NumIntSurf+2, vertexnum, size(case_num_cont, 1));
            T1_RI_SubCortical_temp      = zeros(NumSubSurf, vertexnum, size(case_num_cont, 1));
            T1_PG_IntCortical_temp      = zeros(NumIntSurf+2, vertexnum, size(case_num_cont, 1));
            T1_PG_gw_IntCortical_temp   = zeros(size(case_num_cont, 1), vertexnum);
            T1_PG_SubCortical_temp      = zeros(NumSubSurf, vertexnum, size(case_num_cont, 1));
            T1_TG_IntCortical_temp      = zeros(NumIntSurf+2, vertexnum, size(case_num_cont, 1));
            T1_TG_SubCortical_temp      = zeros(NumSubSurf, vertexnum, size(case_num_cont, 1));
            
            %% Morphological features
            T1_CT_midCortical_temp      = zeros(size(case_num_cont, 1), vertexnum);
            T1_SD_wmCortical_temp       = zeros(size(case_num_cont, 1), vertexnum);
            T1_MC_midCortical_temp      = zeros(size(case_num_cont, 1), vertexnum);
            
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
            file_postfix = [ '_rsl.txt' ];
            for j = 1 : size(case_num_cont, 1)
                fprintf([ case_num_cont{j} ' : ' ]);
                prefix_path = [ data_dir '/' case_num_cont{j} '/measurement/' Prefix_cont '_' case_num_cont{j} ];
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
                    
                    T1_RI_IntCortical_temp(i, 1:vertexnum, j) = SurfStatReadData( { [ basename{1} 'RI_corrected' file_postfix ], [ basename{2} 'RI_corrected' file_postfix ] } );
                    T1_PG_IntCortical_temp(i, 1:vertexnum, j) = SurfStatReadData( { [ basename{1} 'pg' file_postfix ], [ basename{2} 'pg' file_postfix ] } );
                    T1_TG_IntCortical_temp(i, 1:vertexnum, j) = SurfStatReadData( { [ basename{1} 'tg' file_postfix ], [ basename{2} 'tg' file_postfix ] } );
                end
                fprintf([ 'T1: RI_corrected, PG, TG in Inc, ' ]);
                
                %% Morphological / Intensity-based features: MC, SD, PG_GW
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
                
                %% Intensity-based features of subcortical surfaces: RI, pg, tg
                for i = 1 : NumSubSurf
                    basename{1}  = [ prefix_path '_white_surface_' num2str(i) '_left_'  num2str(NumMesh) postfix_surf ];
                    basename{2}  = [ prefix_path '_white_surface_' num2str(i) '_right_'  num2str(NumMesh) postfix_surf ];
                    
                    T1_RI_SubCortical_temp(i, 1:vertexnum, j) = SurfStatReadData( { [ basename{1} 'RI' file_postfix ], [ basename{2} 'RI' file_postfix ] } );
                    T1_PG_SubCortical_temp(i, 1:vertexnum, j) = SurfStatReadData( { [ basename{1} 'pg' file_postfix ], [ basename{2} 'pg' file_postfix ] } );
                    T1_TG_SubCortical_temp(i, 1:vertexnum, j) = SurfStatReadData( { [ basename{1} 'tg' file_postfix ], [ basename{2} 'tg' file_postfix ] } );
                end
                fprintf([ 'RI, PG, TG in Sub\n' ]);
            end
            
        end
        
        % FLAIR
        for FLAIR_read_file = 1
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% FLAIR
            %% Intensity-based features
            FLAIR_RI_IntCortical_temp      = zeros(NumIntSurf+2, vertexnum, size(case_num_cont, 1));
            FLAIR_RI_SubCortical_temp      = zeros(NumSubSurf, vertexnum, size(case_num_cont, 1));
            FLAIR_PG_IntCortical_temp      = zeros(NumIntSurf+2, vertexnum, size(case_num_cont, 1));
            FLAIR_PG_gw_IntCortical_temp   = zeros(size(case_num_cont, 1), vertexnum);
            FLAIR_PG_SubCortical_temp      = zeros(NumSubSurf, vertexnum, size(case_num_cont, 1));
            FLAIR_TG_IntCortical_temp      = zeros(NumIntSurf+2, vertexnum, size(case_num_cont, 1));
            FLAIR_TG_SubCortical_temp      = zeros(NumSubSurf, vertexnum, size(case_num_cont, 1));
            
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
            file_postfix = [ '_rsl.txt' ];
            for j = 1 : size(case_num_cont, 1)
                fprintf([ case_num_cont{j} ' : ' ]);
                prefix_path = [ data_dir '/' case_num_cont{j} '/measurement/' Prefix_cont '_' case_num_cont{j} ];
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
        
    end
    
    % pateints
    for patient_computation = 1
        
        % read T1 files
        for T1_read_file = 1
            
            %% Intensity-based features
            T1_RI_IntCortical_temp2      = zeros(NumIntSurf+2, vertexnum, size(case_num_pat, 1));
            T1_RI_SubCortical_temp2      = zeros(NumSubSurf, vertexnum, size(case_num_pat, 1));
            T1_PG_IntCortical_temp2      = zeros(NumIntSurf+2, vertexnum, size(case_num_pat, 1));
            T1_PG_gw_IntCortical_temp2   = zeros(size(case_num_pat, 1), vertexnum);
            T1_PG_SubCortical_temp2      = zeros(NumSubSurf, vertexnum, size(case_num_pat, 1));
            T1_TG_IntCortical_temp2      = zeros(NumIntSurf+2, vertexnum, size(case_num_pat, 1));
            T1_TG_SubCortical_temp2      = zeros(NumSubSurf, vertexnum, size(case_num_pat, 1));
            
            %% Morphological features
            T1_CT_midCortical_temp2      = zeros(size(case_num_pat, 1), vertexnum);
            T1_SD_wmCortical_temp2       = zeros(size(case_num_pat, 1), vertexnum);
            T1_MC_midCortical_temp2      = zeros(size(case_num_pat, 1), vertexnum);
            
            %% FLAIR
            %% Intensity-based features
            FLAIR_RI_IntCortical_temp2      = zeros(NumIntSurf+2, vertexnum, size(case_num_pat, 1));
            FLAIR_RI_SubCortical_temp2      = zeros(NumSubSurf, vertexnum, size(case_num_pat, 1));
            FLAIR_PG_IntCortical_temp2      = zeros(NumIntSurf+2, vertexnum, size(case_num_pat, 1));
            FLAIR_PG_gw_IntCortical_temp2   = zeros(size(case_num_pat, 1), vertexnum);
            FLAIR_PG_SubCortical_temp2      = zeros(NumSubSurf, vertexnum, size(case_num_pat, 1));
            FLAIR_TG_IntCortical_temp2      = zeros(NumIntSurf+2, vertexnum, size(case_num_pat, 1));
            FLAIR_TG_SubCortical_temp2      = zeros(NumSubSurf, vertexnum, size(case_num_pat, 1));
            
            data_dir = '/local_raid/seokjun/01_project/04_IntracorticalAnalysis/03_result/';
            file_postfix = [ '_rsl.txt' ];
            for j = 1 : size(case_num_pat, 1)
                fprintf([ case_num_pat{j} ' : ' ]);
                prefix_path = [ data_dir '/' case_num_pat{j} '/measurement/' Prefix_pat '_' case_num_pat{j} ];
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
            
            for j = 1 : size(case_num_pat, 1)
                fprintf([ case_num_pat{j} ' : ' ]);
                prefix_path = [ data_dir '/' case_num_pat{j} '/measurement/' Prefix_pat '_' case_num_pat{j} ];
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
        for T1_zscore = 1
            
            
            %% z-score
            T1_RI_IntCortical_z         = zeros(NumIntSurf+2, vertexnum, size(case_num_pat, 1));
            T1_RI_SubCortical_z         = zeros(NumSubSurf,   vertexnum, size(case_num_pat, 1));
            T1_PG_IntCortical_z         = zeros(NumIntSurf+2, vertexnum, size(case_num_pat, 1));
            T1_PG_SubCortical_z         = zeros(NumSubSurf,   vertexnum, size(case_num_pat, 1));
            T1_TG_IntCortical_z         = zeros(NumIntSurf+2, vertexnum, size(case_num_pat, 1));
            T1_TG_SubCortical_z         = zeros(NumSubSurf,   vertexnum, size(case_num_pat, 1));
            T1_CT_midCortical_z         = zeros(size(case_num_pat, 1), vertexnum);
            T1_MC_midCortical_z         = zeros(size(case_num_pat, 1), vertexnum);
            T1_SD_wmCortical_z          = zeros(size(case_num_pat, 1), vertexnum);
            T1_PG_gw_IntCortical_z      = zeros(size(case_num_pat, 1), vertexnum);
            
            FLAIR_RI_IntCortical_z         = zeros(NumIntSurf+2, vertexnum, size(case_num_pat, 1));
            FLAIR_RI_SubCortical_z         = zeros(NumSubSurf,   vertexnum, size(case_num_pat, 1));
            FLAIR_PG_IntCortical_z         = zeros(NumIntSurf+2, vertexnum, size(case_num_pat, 1));
            FLAIR_PG_SubCortical_z         = zeros(NumSubSurf,   vertexnum, size(case_num_pat, 1));
            FLAIR_TG_IntCortical_z         = zeros(NumIntSurf+2, vertexnum, size(case_num_pat, 1));
            FLAIR_TG_SubCortical_z         = zeros(NumSubSurf,   vertexnum, size(case_num_pat, 1));
            FLAIR_PG_gw_IntCortical_z      = zeros(size(case_num_pat, 1), vertexnum);
            
            %% T1: calculate z-score of patients w.r.t control's distribution
            for j = 1 : size(case_num_pat, 1)
                fprintf([ case_num_pat{j} ' : ' ]);
                
                %% Intensity-based features: corrected RI, pg, tg
                T1_RI_IntCortical_z(:, 1:vertexnum, j) = (T1_RI_IntCortical_temp2(:, 1:vertexnum, j) - T1_RI_IntCortical_mean) ./ T1_RI_IntCortical_std;
                T1_PG_IntCortical_z(:, 1:vertexnum, j) = (T1_PG_IntCortical_temp2(:, 1:vertexnum, j) - T1_PG_IntCortical_mean) ./ T1_PG_IntCortical_std;
                T1_TG_IntCortical_z(:, 1:vertexnum, j) = (T1_TG_IntCortical_temp2(:, 1:vertexnum, j) - T1_TG_IntCortical_mean) ./ T1_TG_IntCortical_std;
                T1_RI_IntCortical_z(isnan(T1_RI_IntCortical_z)) = 0; T1_RI_IntCortical_z(isinf(T1_RI_IntCortical_z)) = 0;
                T1_PG_IntCortical_z(isnan(T1_PG_IntCortical_z)) = 0; T1_PG_IntCortical_z(isinf(T1_PG_IntCortical_z)) = 0;
                T1_TG_IntCortical_z(isnan(T1_TG_IntCortical_z)) = 0; T1_PG_IntCortical_z(isinf(T1_PG_IntCortical_z)) = 0;
                fprintf([ 'z-score of T1: RI_corrected, PG, TG in Inc, ' ]);
                
                %% Morphological / Intensity-based features: MC, SD, PG_GW
                T1_PG_gw_IntCortical_z(j, 1:vertexnum) = (T1_PG_gw_IntCortical_temp2(j, 1:vertexnum) - T1_PG_gw_IntCortical_mean) ./ T1_PG_gw_IntCortical_std;
                T1_CT_midCortical_z(j, 1:vertexnum)    = (T1_CT_midCortical_temp2(j, 1:vertexnum) - T1_CT_midCortical_mean) ./ T1_CT_midCortical_std;
                T1_MC_midCortical_z(j, 1:vertexnum)    = (T1_MC_midCortical_temp2(j, 1:vertexnum) - T1_MC_midCortical_mean) ./ T1_MC_midCortical_std;
                T1_SD_wmCortical_z(j, 1:vertexnum)     = (T1_SD_wmCortical_temp2(j, 1:vertexnum) - T1_SD_wmCortical_mean) ./ T1_SD_wmCortical_std;
                T1_PG_gw_IntCortical_z(isnan(T1_PG_gw_IntCortical_z)) = 0; T1_PG_gw_IntCortical_z(isinf(T1_PG_gw_IntCortical_z)) = 0;
                T1_CT_midCortical_z(isnan(T1_CT_midCortical_z)) = 0; T1_CT_midCortical_z(isinf(T1_CT_midCortical_z)) = 0;
                T1_MC_midCortical_z(isnan(T1_MC_midCortical_z)) = 0; T1_MC_midCortical_z(isinf(T1_MC_midCortical_z)) = 0;
                T1_SD_wmCortical_z(isnan(T1_SD_wmCortical_z)) = 0; T1_SD_wmCortical_z(isinf(T1_SD_wmCortical_z)) = 0;
                fprintf([ 'CT, MC, SD, PG_gw in Inc, ' ]);
                
                %% Intensity-based features of subcortical surfaces: RI, pg, tg
                T1_RI_SubCortical_z(:, 1:vertexnum, j) = (T1_RI_SubCortical_temp2(:, 1:vertexnum, j) - T1_RI_SubCortical_mean) ./ T1_RI_SubCortical_std;
                T1_PG_SubCortical_z(:, 1:vertexnum, j) = (T1_PG_SubCortical_temp2(:, 1:vertexnum, j) - T1_PG_SubCortical_mean) ./ T1_PG_SubCortical_std;
                T1_TG_SubCortical_z(:, 1:vertexnum, j) = (T1_TG_SubCortical_temp2(:, 1:vertexnum, j) - T1_TG_SubCortical_mean) ./ T1_TG_SubCortical_std;
                T1_RI_SubCortical_z(isnan(T1_RI_SubCortical_z)) = 0; T1_RI_SubCortical_z(isinf(T1_RI_SubCortical_z)) = 0;
                T1_PG_SubCortical_z(isnan(T1_PG_SubCortical_z)) = 0; T1_PG_SubCortical_z(isinf(T1_PG_SubCortical_z)) = 0;
                T1_TG_SubCortical_z(isnan(T1_TG_SubCortical_z)) = 0; T1_TG_SubCortical_z(isinf(T1_TG_SubCortical_z)) = 0;
                fprintf([ 'RI, PG, TG in Sub\n' ]);
            end
            
        end
        
        % compute FLAIR z-score
        for FLAIR_zscore = 1
            
            %% FLAIR: calculate z-score of patients w.r.t control's distribution
            for j = 1 : size(case_num_pat, 1)
                fprintf([ case_num_pat{j} ' : ' ]);
                
                %% Intensity-based features: corrected RI, pg, tg
                FLAIR_RI_IntCortical_z(:, 1:vertexnum, j) = (FLAIR_RI_IntCortical_temp2(:, 1:vertexnum, j) - FLAIR_RI_IntCortical_mean) ./ FLAIR_RI_IntCortical_std;
                FLAIR_PG_IntCortical_z(:, 1:vertexnum, j) = (FLAIR_PG_IntCortical_temp2(:, 1:vertexnum, j) - FLAIR_PG_IntCortical_mean) ./ FLAIR_PG_IntCortical_std;
                FLAIR_TG_IntCortical_z(:, 1:vertexnum, j) = (FLAIR_TG_IntCortical_temp2(:, 1:vertexnum, j) - FLAIR_TG_IntCortical_mean) ./ FLAIR_TG_IntCortical_std;
                FLAIR_RI_IntCortical_z(isnan(FLAIR_RI_IntCortical_z)) = 0; FLAIR_RI_IntCortical_z(isinf(FLAIR_RI_IntCortical_z)) = 0;
                FLAIR_PG_IntCortical_z(isnan(FLAIR_PG_IntCortical_z)) = 0; FLAIR_PG_IntCortical_z(isinf(FLAIR_PG_IntCortical_z)) = 0;
                FLAIR_TG_IntCortical_z(isnan(FLAIR_TG_IntCortical_z)) = 0; FLAIR_TG_IntCortical_z(isinf(FLAIR_TG_IntCortical_z)) = 0;
                fprintf([ 'z-score of FLAIR: RI_corrected, PG, TG in Inc, ' ]);
                
                %% Morphological / Intensity-based features: MC, SD, PG_GW
                FLAIR_PG_gw_IntCortical_z(j, 1:vertexnum) = (FLAIR_PG_gw_IntCortical_temp2(j, 1:vertexnum) - FLAIR_PG_gw_IntCortical_mean) ./ FLAIR_PG_gw_IntCortical_std;
                FLAIR_PG_gw_IntCortical_z(isnan(FLAIR_PG_gw_IntCortical_z)) = 0;
                fprintf([ 'PG_gw in Inc, ' ]);
                
                %% Intensity-based features of subcortical surfaces: RI, pg, tg
                FLAIR_RI_SubCortical_z(:, 1:vertexnum, j) = (FLAIR_RI_SubCortical_temp2(:, 1:vertexnum, j) - FLAIR_RI_SubCortical_mean) ./ FLAIR_RI_SubCortical_std;
                FLAIR_PG_SubCortical_z(:, 1:vertexnum, j) = (FLAIR_PG_SubCortical_temp2(:, 1:vertexnum, j) - FLAIR_PG_SubCortical_mean) ./ FLAIR_PG_SubCortical_std;
                FLAIR_TG_SubCortical_z(:, 1:vertexnum, j) = (FLAIR_TG_SubCortical_temp2(:, 1:vertexnum, j) - FLAIR_TG_SubCortical_mean) ./ FLAIR_TG_SubCortical_std;
                FLAIR_RI_SubCortical_z(isnan(FLAIR_RI_SubCortical_z)) = 0; FLAIR_RI_SubCortical_z(isinf(FLAIR_RI_SubCortical_z)) = 0;
                FLAIR_PG_SubCortical_z(isnan(FLAIR_PG_SubCortical_z)) = 0; FLAIR_PG_SubCortical_z(isinf(FLAIR_PG_SubCortical_z)) = 0;
                FLAIR_TG_SubCortical_z(isnan(FLAIR_TG_SubCortical_z)) = 0; FLAIR_TG_SubCortical_z(isinf(FLAIR_TG_SubCortical_z)) = 0;
                fprintf([ 'RI, PG, TG in Sub\n' ]);
            end
            
        end
   
        % save zscore as a mat file
        for save_zscore = 1
            if(~exist([ 'zscore_database_perilesional_T1_FLAIR_nosmoothing.mat'], 'file'))
                
                save('zscore_database_perilesional_T1_FLAIR_nosmoothing.mat', ...
                    'T1_RI_IntCortical_z', ...
                    'T1_PG_IntCortical_z', ...
                    'T1_PG_gw_IntCortical_z', ...
                    'T1_TG_IntCortical_z', ...
                    'T1_RI_SubCortical_z', ...
                    'T1_CT_midCortical_z', ...
                    'T1_MC_midCortical_z', ...
                    'T1_SD_wmCortical_z', ...
                    'FLAIR_RI_IntCortical_z', ...
                    'FLAIR_PG_IntCortical_z', ...
                    'FLAIR_PG_gw_IntCortical_z', ...
                    'FLAIR_TG_IntCortical_z', ...
                    'FLAIR_RI_SubCortical_z');
                
            end
            
        end

    end

    % surface visualization mapping MRI values
    for visualization = 1
        
        %% Geodesic distance for the lesion
        OUTPATH          = '/local_raid/seokjun/01_project/04_IntracorticalAnalysis/03_result/';
        lesion_label_dir = '/local_raid/CTFCD-1.2.0_64/Lesion/lesion_surf/';
        postfix_surf     = '_native';
        geoDist    = cell(size(case_num_pat, 1), 1);
        surf_ind   = cell(size(case_num_pat, 1), 1);
        surf   = SurfStatReadSurf({ '/local_raid/seokjun/01_project/04_IntracorticalAnalysis/03_result/002/surfaces/mcd_002_gray_surface_left_81920_native_t1.obj', ...
                                    '/local_raid/seokjun/01_project/04_IntracorticalAnalysis/03_result/002/surfaces/mcd_002_gray_surface_right_81920_native_t1.obj' });
        surf = surfGetNeighborsHong(surf);
        nbr = surf.nbr';
        
        for j = 1 : size(case_num_pat, 1)
            idx = find(strcmp(case_num_pat, case_num_pat{j}));
            
            if(exist([ lesion_label_dir '/' Prefix_pat '_' case_num_pat{j} '_label_union_left_rsl.txt' ], 'file'))
                lesion_label_data = SurfStatReadData( { [ lesion_label_dir '/' Prefix_pat '_' case_num_pat{j} '_label_union_left_rsl.txt' ], ...
                    [ lesion_label_dir '/' Prefix_pat '_' case_num_pat{j} '_label_union_right_rsl.txt' ] } );
                
                geoDist_intra{j}.surf = zeros(NumIntSurf + 2 + NumSubSurf, 81924);
                surf_ind{j}.surf = cell(NumIntSurf + 2 + NumSubSurf, 1);
                for i = 1 : NumIntSurf + 2
                    if(i == 1)
                        basename  = [ Prefix_pat '_' case_num_pat{j} '_gray_surface' ];
                        temp_str = 'gm';
                    elseif(i == NumIntSurf + 2)
                        basename  = [ Prefix_pat '_' case_num_pat{j} '_white_surface' ];
                        temp_str = 'wm';
                    else
                        basename  = [ Prefix_pat '_' case_num_pat{j} '_intracortical_surface_' num2str(i-1) ];
                        temp_str = num2str(i-1);
                    end
                    
                    surf_l_name = [ OUTPATH case_num_pat{j} '/surfaces/' basename '_left_'  num2str(NumMesh) postfix_surf '_t1.obj' ];
                    surf_r_name = [ OUTPATH case_num_pat{j} '/surfaces/' basename '_right_' num2str(NumMesh) postfix_surf '_t1.obj' ];
                    surf_l_rsl_name = [ OUTPATH case_num_pat{j} '/surfaces/' basename '_left_'  num2str(NumMesh) postfix_surf '_t1_rsl.obj' ];
                    surf_r_rsl_name = [ OUTPATH case_num_pat{j} '/surfaces/' basename '_right_' num2str(NumMesh) postfix_surf '_t1_rsl.obj' ];
                    xfm_l_name  = [ OUTPATH case_num_pat{j} '/xfm/' Prefix_pat '_' case_num_pat{j} '_left_surfmap.sm' ];
                    xfm_r_name  = [ OUTPATH case_num_pat{j} '/xfm/' Prefix_pat '_' case_num_pat{j} '_right_surfmap.sm' ];
                    
                    if(exist(surf_l_rsl_name, 'file') & exist(surf_r_rsl_name, 'file'))
                        surf   = SurfStatReadSurf({ surf_l_rsl_name, surf_l_rsl_name });
                        surf_ind{j}.surf{i}.tri = surf.tri;
                        surf_ind{j}.surf{i}.coord = surf.coord;
                        surf.nbr = nbr;
                        geoDist{j}.surf(i, :) = surfGeoDistHong(surf, lesion_label_data);
                        fprintf('%s ', temp_str);
                        continue;
                    end
                    
                    if(~exist(surf_l_name, 'file') | ~exist(surf_r_name, 'file'))
                        continue;
                    end
                    
                    [r, s] = system(['sphere_resample_obj ' surf_l_name ' ' xfm_l_name ' ' surf_l_rsl_name ]);
                    [r, s] = system(['sphere_resample_obj ' surf_r_name ' ' xfm_r_name ' ' surf_r_rsl_name ]);
                    
                    surf   = SurfStatReadSurf({ surf_l_rsl_name, surf_l_rsl_name });
                    surf_ind{j}.surf{i}.tri = surf.tri;
                    surf_ind{j}.surf{i}.coord = surf.coord;
                    surf.nbr = nbr;
                    geoDist{j}.surf(i, :) = surfGeoDistHong(surf, lesion_label_data);
                    fprintf('%s ', temp_str);
                end
                for i = 1 : NumSubSurf
                    basename = [ Prefix_pat '_' case_num_pat{j} '_white_surface_' num2str(i) ];
                    temp_str = num2str(i);
                    
                    surf_l_name = [ OUTPATH case_num_pat{j} '/surfaces/' basename '_left_'  num2str(NumMesh) postfix_surf '_t1.obj' ];
                    surf_r_name = [ OUTPATH case_num_pat{j} '/surfaces/' basename '_right_' num2str(NumMesh) postfix_surf '_t1.obj' ];
                    surf_l_rsl_name = [ OUTPATH case_num_pat{j} '/surfaces/' basename '_left_'  num2str(NumMesh) postfix_surf '_t1_rsl.obj' ];
                    surf_r_rsl_name = [ OUTPATH case_num_pat{j} '/surfaces/' basename '_right_' num2str(NumMesh) postfix_surf '_t1_rsl.obj' ];
                    xfm_l_name  = [ OUTPATH case_num_pat{j} '/xfm/' Prefix_pat '_' case_num_pat{j} '_left_surfmap.sm' ];
                    xfm_r_name  = [ OUTPATH case_num_pat{j} '/xfm/' Prefix_pat '_' case_num_pat{j} '_right_surfmap.sm' ];
                    
                    if(exist(surf_l_rsl_name, 'file') & exist(surf_r_rsl_name, 'file'))
                        surf   = SurfStatReadSurf({ surf_l_rsl_name, surf_l_rsl_name });
                        surf_ind{j}.surf{i+5}.tri = surf.tri;
                        surf_ind{j}.surf{i+5}.coord = surf.coord;
                        
                        surf.nbr = nbr;
                        geoDist{j}.surf(i+5, :) = surfGeoDistHong(surf, lesion_label_data);
                        fprintf('%s ', temp_str);
                        continue;
                    end
                    
                    if(~exist(surf_l_name, 'file') | ~exist(surf_r_name, 'file'))
                        continue;
                    end
                    
                    [r, s] = system(['sphere_resample_obj ' surf_l_name ' ' xfm_l_name ' ' surf_l_rsl_name ]);
                    [r, s] = system(['sphere_resample_obj ' surf_r_name ' ' xfm_r_name ' ' surf_r_rsl_name ]);
                    
                    surf   = SurfStatReadSurf({ surf_l_rsl_name, surf_l_rsl_name });
                    surf_ind{j}.surf{i+5}.tri = surf.tri;
                    surf_ind{j}.surf{i+5}.coord = surf.coord;
                    
                    surf.nbr = nbr;
                    geoDist{j}.surf(i+5, :) = surfGeoDistHong(surf, lesion_label_data);
                    fprintf('%s ', temp_str);
                end
            end
            disp(['case : ' case_num_pat{j} ' done']);
        end
        
    end
    
end
