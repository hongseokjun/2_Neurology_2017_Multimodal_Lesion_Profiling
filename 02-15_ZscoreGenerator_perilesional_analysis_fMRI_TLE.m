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
    NumSubSurf                = 6
    NumMesh                   = 81920
    SamplingSpace             = 'native'
    Kernel                    = 2
    weight_degree             = 1
    Parametric                = 'quadratic'
    average_surface_dir       = '/local_raid/seokjun/01_project/04_IntracorticalAnalysis/03_result/average_surfaces/'
    
    visualization = 1;
    load('/local_raid/seokjun/01_project/04_IntracorticalAnalysis/01_analysis/pipeline/postprocessing_zscore_visualization/colormap_noel1.mat');
    
    fid = fopen(Cases_healthy_cont);
    demo = textscan(fid, '%s%d%s%s%s', 'Delimiter', ',', 'CollectOutput', 1);
    case_num_healthy_cont = demo{1};
    age_healthy_cont = demo{2};
    gender_healthy_cont = demo{3}(:, 1);
    fclose(fid);
    
    fid = fopen(Cases_disease_cont);
    demo = textscan(fid, '%s%f%s%s', 'Delimiter', ',', 'CollectOutput', 1);
    case_num_disease_cont       = demo{1};
    age_disease_cont            = demo{2};
    gender_disease_cont         = demo{3}(:, 1);
    seizure_lateralization  = demo{3}(:, 2);
    fclose(fid);    
    
    fid = fopen(Cases_pat);
    demo = textscan(fid, '%s%d%s%s%s', 'Delimiter', ',', 'CollectOutput', 1);
    case_num_pat = demo{1};
    age_pat      = demo{2};
    gender_pat   = demo{3}(:, 1);
    histo_pat    = demo{3}(:, 3);
    
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
    
    % Modality: resting-fMRI
    
    %% Feature: ALFF (numControl x vertexnum),
    %%          ReHo (numControl x vertexnum),
    
    %% Abb: ALFF Amplitude of low frequency fluctuation
    %%      ReHo Regional homogeneity

    % healthy controls
    for healthy_control_computation = 1
        
        % fMRI
        for fMRI_read_file = 1
            
            fMRI_ALFF_temp = zeros(size(case_num_healthy_cont, 1), vertexnum);
            fMRI_ReHo_temp = zeros(size(case_num_healthy_cont, 1), vertexnum);
            
            data_dir = '/local_raid/seokjun/01_project/05_NoelRest/03_Result/01_DPARSF_processing/DPARSFA_scrub/02_without_mean_sig_regout/bbr_rest/';
            file_postfix = [ '_aniso_sm_' num2str(Kernel) '_rsl.txt' ];
            for j = 1 : size(case_num_healthy_cont, 1)
                fprintf([ case_num_healthy_cont{j} ' : ' ]);
                prefix_path = [ data_dir '/' Prefix_healthy_cont '_' case_num_healthy_cont{j} '/' Prefix_healthy_cont '_' case_num_healthy_cont{j} '_z' ];                
                                     
                fMRI_ALFF_temp(j, 1:vertexnum) = SurfStatReadData( { [ prefix_path 'ALFFMap_mid_left_' num2str(NumMesh) file_postfix ], [ prefix_path 'ALFFMap_mid_right_' num2str(NumMesh) file_postfix ] } );
                fMRI_ReHo_temp(j, 1:vertexnum) = SurfStatReadData( { [ prefix_path 'ReHoMap_mid_left_' num2str(NumMesh) file_postfix ], [ prefix_path 'ReHoMap_mid_right_' num2str(NumMesh) file_postfix ] } );
                fprintf([ 'rs-fMRI: ALFF and ReHo\n' ]);
                
            end
            
        end
    
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
            fMRI_ALFF_mean = mean(fMRI_ALFF_temp(:, 1:vertexnum), 1);
            fMRI_ReHo_mean = mean(fMRI_ReHo_temp(:, 1:vertexnum), 1);
            fMRI_ALFF_std  = std(fMRI_ALFF_temp(:, 1:vertexnum), 0, 1);
            fMRI_ReHo_std  = std(fMRI_ReHo_temp(:, 1:vertexnum), 0, 1);
            
        end
        
                % compute fMRI z-score
        for fMRI_zscore = 1
            
            %% z-score
            fMRI_ALFF_z_cont = zeros(size(case_num_healthy_cont, 1), vertexnum);
            fMRI_ReHo_z_cont = zeros(size(case_num_healthy_cont, 1), vertexnum);
            
            %% fMRI: calculate z-score of patients w.r.t control's distribution
            for j = 1 : size(case_num_healthy_cont, 1)
                fprintf([ case_num_healthy_cont{j} ' : ' ]);
                
                %% Intensity-based features: corrected RI, pg, tg
                fMRI_ALFF_z_cont(j, 1:vertexnum) = (fMRI_ALFF_temp(j, 1:vertexnum) - fMRI_ALFF_mean) ./ fMRI_ALFF_std;
                fMRI_ReHo_z_cont(j, 1:vertexnum) = (fMRI_ReHo_temp(j, 1:vertexnum) - fMRI_ReHo_mean) ./ fMRI_ReHo_std;
                fprintf([ 'z-score of fMRI: ALFF, ReHo\n' ]);
            end
            
            fMRI_ALFF_z_cont(find(isinf(fMRI_ALFF_z_cont))) = 0;
            fMRI_ReHo_z_cont(find(isinf(fMRI_ReHo_z_cont))) = 0;
            
%             fMRI_ALFF_z_cont(logical(outlier_vector), :) = [];
%             fMRI_ReHo_z_cont(logical(outlier_vector), :) = [];
            
        end
    end
    
    % disease controls
    for disease_control_computation = 1
        
        % fMRI
        for fMRI_read_file = 1
            
            fMRI_ALFF_temp2 = zeros(size(case_num_disease_cont, 1), vertexnum);
            fMRI_ReHo_temp2 = zeros(size(case_num_disease_cont, 1), vertexnum);
            
            data_dir = '/host/weka/export02/data/min/fMRI/DPARSF_processing/DPARSFA_scrub/02_without_mean_sig_regout/bbr_rest/';
            file_postfix = [ '_aniso_sm_' num2str(Kernel) '_rsl.txt' ];
            for j = 1 : size(case_num_disease_cont, 1)
                fprintf([ case_num_disease_cont{j} ' : ' ]);
                prefix_path = [ data_dir '/' Prefix_disease_cont '_' case_num_disease_cont{j} '/' Prefix_disease_cont '_' case_num_disease_cont{j} '_z' ];                
                                     
                fMRI_ALFF_temp2(j, 1:vertexnum) = SurfStatReadData( { [ prefix_path 'ALFFMap_mid_left_' num2str(NumMesh) file_postfix ], [ prefix_path 'ALFFMap_mid_right_' num2str(NumMesh) file_postfix ] } );
                fMRI_ReHo_temp2(j, 1:vertexnum) = SurfStatReadData( { [ prefix_path 'ReHoMap_mid_left_' num2str(NumMesh) file_postfix ], [ prefix_path 'ReHoMap_mid_right_' num2str(NumMesh) file_postfix ] } );
                fprintf([ 'rs-fMRI: ALFF and ReHo\n' ]);
                
            end
            
        end            
                % compute fMRI z-score
        for fMRI_zscore = 1
            
            %% z-score
            fMRI_ALFF_z_disease_cont = zeros(size(case_num_disease_cont, 1), vertexnum);
            fMRI_ReHo_z_disease_cont = zeros(size(case_num_disease_cont, 1), vertexnum);
            
            %% fMRI: calculate z-score of patients w.r.t control's distribution
            for j = 1 : size(case_num_disease_cont, 1)
                fprintf([ case_num_healthy_cont{j} ' : ' ]);
                
                %% Intensity-based features: corrected RI, pg, tg
                fMRI_ALFF_z_disease_cont(j, 1:vertexnum) = (fMRI_ALFF_temp2(j, 1:vertexnum) - fMRI_ALFF_mean) ./ fMRI_ALFF_std;
                fMRI_ReHo_z_disease_cont(j, 1:vertexnum) = (fMRI_ReHo_temp2(j, 1:vertexnum) - fMRI_ReHo_mean) ./ fMRI_ReHo_std;
                fprintf([ 'z-score of fMRI: ALFF, ReHo\n' ]);
            end
            
            fMRI_ALFF_z_disease_cont(find(isinf(fMRI_ALFF_z_disease_cont))) = 0;
            fMRI_ReHo_z_disease_cont(find(isinf(fMRI_ReHo_z_disease_cont))) = 0;

        end
    end
    
end
