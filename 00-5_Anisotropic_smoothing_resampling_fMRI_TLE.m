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
    NumSubSurf                = 3
    NumMesh                   = 81920
    SamplingSpace             = 'native'
    img_contrast              = {'t1'}
    Kernel                    = 2
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
    
    LESIONPATH = '/data/noel/noel6/CTFCD-1.2.0_64/Lesion/lesion_surf/';    
    
    %% Features: RI, RI_corrected, pg, pg_GM_WM, tg, ct, sd2, mc, FA, MD
    %% Modalities: T1 FLAIR DTI, T1_FLAIR_ratio
    %% Surfaces: intra 1,2,3,4,5 | subcortical 1,2,3,4,5,6
    Feature    = { 'ALFF', 'ReHo' };
    Modality   = { 'rest' };
    
    if(strcmp(SamplingSpace, 'native'))
        space='native';
        postfix_surf='_native';
        postfix_MRI='_nuc';
    elseif(strcmp(SamplingSpace, 'tal'))
        space='final';
        postfix_surf='';
        postfix_MRI='_final';
    end
    
    %% feature-specific parameter tuning (weight degree for extralesion and lesion)
    aniso_smooth_parameters = { [ 1.5 3 ],          %% ALFF
                                [ 1.5 3 ] };        %% ReHo

end

%% Smooth the data outside the given mask to be unaffected by the abnormalities in the lesion
for anisotropic_smooth = 1
    
    % Healthy control
    % In controls, as there is no lesion mask, we just smoothed whole
    % vertices using anisotropic kernel which gives almost same result of
    % heat kernel smoothing at this time.
    OUTPATH = '/local_raid/seokjun/01_project/05_NoelRest/03_Result/01_DPARSF_processing/DPARSFA_scrub/02_without_mean_sig_regout/bbr_rest/';
    for k = 1 : size(case_num_healthy_cont, 1)
        fprintf('Control case: %s start! ', case_num_healthy_cont{k});
        
        tic
        [r, s] = system(['gunzip -vf ' OUTPATH '/' Prefix_healthy_cont '_' case_num_healthy_cont{k} '/' Prefix_healthy_cont '_' case_num_healthy_cont{k} '*.obj.gz']);
        surf_l_name = [ OUTPATH '/' Prefix_healthy_cont '_' case_num_healthy_cont{k} '/' Prefix_healthy_cont '_' case_num_healthy_cont{k} '_mid_surface_left_' num2str(NumMesh) '.obj' ];
        surf_r_name = [ OUTPATH '/' Prefix_healthy_cont '_' case_num_healthy_cont{k} '/' Prefix_healthy_cont '_' case_num_healthy_cont{k} '_mid_surface_right_' num2str(NumMesh) '.obj' ];
        surf_final = SurfStatReadSurf({ surf_l_name, surf_r_name });
        MidMask     = SurfStatMaskCut(surf_final);
        for m = 1 : size(Modality, 2)
            Modality_temp = Modality{m};
           
            surf_l_name = [ OUTPATH '/' Prefix_healthy_cont '_' case_num_healthy_cont{k} '/'  Prefix_healthy_cont '_' case_num_healthy_cont{k} '_mid_surface_left_' num2str(NumMesh) '_' Modality_temp '_bbr.obj' ];
            surf_r_name = [ OUTPATH '/' Prefix_healthy_cont '_' case_num_healthy_cont{k} '/' Prefix_healthy_cont '_' case_num_healthy_cont{k} '_mid_surface_right_' num2str(NumMesh) '_' Modality_temp '_bbr.obj' ];
            
            surf_l = SurfStatReadSurf(surf_l_name);
            surf_r = SurfStatReadSurf(surf_r_name);
            surf   = SurfStatReadSurf({ surf_l_name, surf_r_name });
                
            for f = 1 : size(Feature, 2)
                Y_l_name = [ OUTPATH '/' Prefix_healthy_cont '_' case_num_healthy_cont{k} '/'  Prefix_healthy_cont '_' case_num_healthy_cont{k} '_z'  Feature{f} 'Map_mid_left_' num2str(NumMesh) '.txt' ];
                Y_r_name = [ OUTPATH '/' Prefix_healthy_cont '_' case_num_healthy_cont{k} '/'  Prefix_healthy_cont '_' case_num_healthy_cont{k} '_z'  Feature{f} 'Map_mid_right_' num2str(NumMesh) '.txt' ];
                if(exist(Y_l_name, 'file') & exist(Y_r_name, 'file'))
                    weight_degree = aniso_smooth_parameters{f}(1);
                    Y_l = SurfStatReadData(Y_l_name);
                    Y_r = SurfStatReadData(Y_r_name);
                    Y_l(MidMask(1:40962)==0) = 0; Y_r(MidMask(40963:81924)==0) = 0;
                    Y_l_smooth = SurfStatSmoothOutsideMask( Y_l, surf_l, Kernel, weight_degree );
                    Y_r_smooth = SurfStatSmoothOutsideMask( Y_r, surf_r, Kernel, weight_degree );
                    Y_l_smooth_name = [  OUTPATH '/' Prefix_healthy_cont '_' case_num_healthy_cont{k} '/'  Prefix_healthy_cont '_' case_num_healthy_cont{k} '_z'  Feature{f} 'Map_mid_left_' num2str(NumMesh) '_aniso_sm_' num2str(Kernel) '.txt' ];
                    Y_r_smooth_name = [  OUTPATH '/' Prefix_healthy_cont '_' case_num_healthy_cont{k} '/'  Prefix_healthy_cont '_' case_num_healthy_cont{k} '_z'  Feature{f} 'Map_mid_right_' num2str(NumMesh) '_aniso_sm_' num2str(Kernel) '.txt' ];
                    SurfStatWriteData(Y_l_smooth_name, Y_l_smooth);
                    SurfStatWriteData(Y_r_smooth_name, Y_r_smooth);
                end
            end
        end
        
        toc
    end
    
    % Disease control
    % The areas outside and inside the mask are separately smoothed and
    % merged together at the end.
    OUTPATH = '/host/weka/export02/data/min/fMRI/DPARSF_processing/DPARSFA_scrub/02_without_mean_sig_regout/bbr_rest/';
    for k = 1 : size(case_num_disease_cont, 1)
        fprintf('Disease control case: %s start! ', case_num_disease_cont{k});
        
        tic
        [r, s] = system(['gunzip -vf ' OUTPATH '/' Prefix_disease_cont '_' case_num_disease_cont{k} '/' Prefix_disease_cont '_' case_num_disease_cont{k} '*.obj.gz']);
        surf_l_name = [ OUTPATH '/' Prefix_disease_cont '_' case_num_disease_cont{k} '/' Prefix_disease_cont '_' case_num_disease_cont{k} '_mid_surface_left_' num2str(NumMesh) '.obj' ];
        surf_r_name = [ OUTPATH '/' Prefix_disease_cont '_' case_num_disease_cont{k} '/' Prefix_disease_cont '_' case_num_disease_cont{k} '_mid_surface_right_' num2str(NumMesh) '.obj' ];
        surf_final = SurfStatReadSurf({ surf_l_name, surf_r_name });
        MidMask     = SurfStatMaskCut(surf_final);
        for m = 1 : size(Modality, 2)
            Modality_temp = Modality{m};
           
            surf_l_name = [ OUTPATH '/' Prefix_disease_cont '_' case_num_disease_cont{k} '/'  Prefix_disease_cont '_' case_num_disease_cont{k} '_mid_surface_left_' num2str(NumMesh) '_' Modality_temp '_bbr.obj' ];
            surf_r_name = [ OUTPATH '/' Prefix_disease_cont '_' case_num_disease_cont{k} '/' Prefix_disease_cont '_' case_num_disease_cont{k} '_mid_surface_right_' num2str(NumMesh) '_' Modality_temp '_bbr.obj' ];
            
            surf_l = SurfStatReadSurf(surf_l_name);
            surf_r = SurfStatReadSurf(surf_r_name);
            surf   = SurfStatReadSurf({ surf_l_name, surf_r_name });
                
            for f = 1 : size(Feature, 2)
                Y_l_name = [ OUTPATH '/' Prefix_disease_cont '_' case_num_disease_cont{k} '/'  Prefix_disease_cont '_' case_num_disease_cont{k} '_z'  Feature{f} 'Map_mid_left_' num2str(NumMesh) '.txt' ];
                Y_r_name = [ OUTPATH '/' Prefix_disease_cont '_' case_num_disease_cont{k} '/'  Prefix_disease_cont '_' case_num_disease_cont{k} '_z'  Feature{f} 'Map_mid_right_' num2str(NumMesh) '.txt' ];
                if(exist(Y_l_name, 'file') & exist(Y_r_name, 'file'))
                    weight_degree = aniso_smooth_parameters{f}(1);
                    Y_l = SurfStatReadData(Y_l_name);
                    Y_r = SurfStatReadData(Y_r_name);
                    Y_l(MidMask(1:40962)==0) = 0; Y_r(MidMask(40963:81924)==0) = 0;
                    Y_l_smooth = SurfStatSmoothOutsideMask( Y_l, surf_l, Kernel, weight_degree );
                    Y_r_smooth = SurfStatSmoothOutsideMask( Y_r, surf_r, Kernel, weight_degree );
                    Y_l_smooth_name = [  OUTPATH '/' Prefix_disease_cont '_' case_num_disease_cont{k} '/'  Prefix_disease_cont '_' case_num_disease_cont{k} '_z'  Feature{f} 'Map_mid_left_' num2str(NumMesh) '_aniso_sm_' num2str(Kernel) '.txt' ];
                    Y_r_smooth_name = [  OUTPATH '/' Prefix_disease_cont '_' case_num_disease_cont{k} '/'  Prefix_disease_cont '_' case_num_disease_cont{k} '_z'  Feature{f} 'Map_mid_right_' num2str(NumMesh) '_aniso_sm_' num2str(Kernel) '.txt' ];
                    SurfStatWriteData(Y_l_smooth_name, Y_l_smooth);
                    SurfStatWriteData(Y_r_smooth_name, Y_r_smooth);
                end
            end
        end
        
        toc
    end
    
end

%% Resample the features
for resample_features = 1
        
    %% Healthy control resampling ...
    OUTPATH = '/local_raid/seokjun/01_project/05_NoelRest/03_Result/01_DPARSF_processing/DPARSFA_scrub/02_without_mean_sig_regout/bbr_rest/';
    XFMPATH = '/local_raid/seokjun/01_project/04_IntracorticalAnalysis/03_result/';
    for k = 1 : size(case_num_healthy_cont, 1)
        fprintf('Control case: %s start! ', case_num_healthy_cont{k});
        tic
        
        for m = 1 : size(Modality, 2)
            Modality_temp = Modality{m};
                 
            surf_l_xfm = [ XFMPATH case_num_healthy_cont{k} '/xfm/' Prefix_healthy_cont '_' case_num_healthy_cont{k} '_left_surfmap.sm' ];
            surf_r_xfm = [ XFMPATH case_num_healthy_cont{k} '/xfm/' Prefix_healthy_cont '_' case_num_healthy_cont{k} '_right_surfmap.sm' ];
            
            for f = 1 : size(Feature, 2)
                Y_l_smooth_name = [  OUTPATH '/' Prefix_healthy_cont '_' case_num_healthy_cont{k} '/'  Prefix_healthy_cont '_' case_num_healthy_cont{k} '_z'  Feature{f} 'Map_mid_left_' num2str(NumMesh) '_aniso_sm_' num2str(Kernel) ];
                Y_r_smooth_name = [  OUTPATH '/' Prefix_healthy_cont '_' case_num_healthy_cont{k} '/'  Prefix_healthy_cont '_' case_num_healthy_cont{k} '_z'  Feature{f} 'Map_mid_right_' num2str(NumMesh) '_aniso_sm_' num2str(Kernel) ];
                if(exist([Y_l_smooth_name '.txt'], 'file') & exist([Y_r_smooth_name '.txt'], 'file'))
                    [r, s] = system(['/data/noel/noel6/Oct-2010/CIVET-1.2.0/progs/surface-resample2 -clobber ' surf_l_xfm ' ' Y_l_smooth_name '.txt ' Y_l_smooth_name '_rsl.txt'])
                    [r, s] = system(['/data/noel/noel6/Oct-2010/CIVET-1.2.0/progs/surface-resample2 -clobber ' surf_r_xfm ' ' Y_r_smooth_name '.txt ' Y_r_smooth_name '_rsl.txt'])
                end
            end

        end
        toc
    end
    
    %% Disease control resampling ...
    OUTPATH = '/host/weka/export02/data/min/fMRI/DPARSF_processing/DPARSFA_scrub/02_without_mean_sig_regout/bbr_rest/';
    XFMPATH = '/local_raid/seokjun/01_project/04_IntracorticalAnalysis/03_result/97_3T_TLE_HS/';
    for k = 1 : size(case_num_disease_cont, 1)
        fprintf('Disease control case: %s start! ', case_num_disease_cont{k});
        tic
        
        for m = 1 : size(Modality, 2)
           Modality_temp = Modality{m};
                        
            surf_l_xfm = [ XFMPATH case_num_disease_cont{k} '/xfm/' Prefix_disease_cont '_' case_num_disease_cont{k} '_left_surfmap.sm' ];
            surf_r_xfm = [ XFMPATH case_num_disease_cont{k} '/xfm/' Prefix_disease_cont '_' case_num_disease_cont{k} '_right_surfmap.sm' ];
                        
            for f = 1 : size(Feature, 2)
                Y_l_smooth_name = [  OUTPATH '/' Prefix_disease_cont '_' case_num_disease_cont{k} '/'  Prefix_disease_cont '_' case_num_disease_cont{k} '_z'  Feature{f} 'Map_mid_left_' num2str(NumMesh) '_aniso_sm_' num2str(Kernel) ];
                Y_r_smooth_name = [  OUTPATH '/' Prefix_disease_cont '_' case_num_disease_cont{k} '/'  Prefix_disease_cont '_' case_num_disease_cont{k} '_z'  Feature{f} 'Map_mid_right_' num2str(NumMesh) '_aniso_sm_' num2str(Kernel) ];
                if(exist([Y_l_smooth_name '.txt'], 'file') & exist([Y_r_smooth_name '.txt'], 'file'))
                    [r, s] = system(['/data/noel/noel6/Oct-2010/CIVET-1.2.0/progs/surface-resample2 -clobber ' surf_l_xfm ' ' Y_l_smooth_name '.txt ' Y_l_smooth_name '_rsl.txt'])
                    [r, s] = system(['/data/noel/noel6/Oct-2010/CIVET-1.2.0/progs/surface-resample2 -clobber ' surf_r_xfm ' ' Y_r_smooth_name '.txt ' Y_r_smooth_name '_rsl.txt'])
                end
            end
 
        end
        toc
    end
    
end
