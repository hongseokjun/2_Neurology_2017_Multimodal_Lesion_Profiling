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
    NumSubSurf                = 6
    NumMesh                   = 81920
    SamplingSpace             = 'native'
    Kernel                    = 5
    weight_degree             = 1
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
    
end

%% Set up some parameters and variables
for setup_param = 1
    
    OUTPATH = '/local_raid/seokjun/01_project/04_IntracorticalAnalysis/03_result/';
    LESIONPATH = '/local_raid/CTFCD-1.2.0_64/Lesion/lesion_surf/';    
    
    %% Features: RI, RI_corrected, pg, pg_GM_WM, tg, ct, sd2, mc, FA, MD
    %% Modalities: T1 FLAIR DTI, T1_FLAIR_ratio
    %% Surfaces: intra 1,2,3,4,5 | subcortical 1,2,3,4,5,6
    Feature    = { 'RI', 'RI_corrected', 'pg', 'pg_GM_WM', 'tg', 'ct', 'sd2', 'mc', 'FA', 'MD' };
    Modality   = { 't1', 'flair', 'dti', 't1_flair_ratio' };
    
    if(strcmp(SamplingSpace, 'native'))
        space='native';
        postfix_surf='_native';
        postfix_MRI='_nuc';
    elseif(strcmp(SamplingSpace, 'tal'))
        space='final';
        postfix_surf='';
        postfix_MRI='_final';
    end
    
end

%% Resample the features
for resample_features = 1
    
    % control
    for k = 1 : size(case_num_cont, 1)
        fprintf('Control case: %s start! ', case_num_cont{k});
        tic
        
        for m = 1 : size(Modality, 2)
            Modality_temp = Modality{m};
            if(strcmp(Modality_temp, 't1_flair_ratio'))
                Modality_temp = 't1';
            end
            
            surf_l_xfm = [ OUTPATH case_num_cont{k} '/xfm/' Prefix_cont '_' case_num_cont{k} '_left_surfmap.sm' ];
            surf_r_xfm = [ OUTPATH case_num_cont{k} '/xfm/' Prefix_cont '_' case_num_cont{k} '_right_surfmap.sm' ];
            
            for i = 1 : NumIntSurf + 2
                if(i == 1)
                    basename  = [ Prefix_cont '_' case_num_cont{k} '_gray_surface' ];
                elseif(i == NumIntSurf + 2)
                    basename  = [ Prefix_cont '_' case_num_cont{k} '_white_surface' ];
                else
                    basename  = [ Prefix_cont '_' case_num_cont{k} '_intracortical_surface_' num2str(i-1) ];
                end
                
                for f = 1 : size(Feature, 2)
                    Y_l_smooth_name = [ OUTPATH case_num_cont{k} '/measurement/' basename '_left_'  num2str(NumMesh) postfix_surf '_' Modality{m} '_' Feature{f} ];
                    Y_r_smooth_name = [ OUTPATH case_num_cont{k} '/measurement/' basename '_right_'  num2str(NumMesh) postfix_surf '_' Modality{m} '_' Feature{f} ];
                    if(exist([Y_l_smooth_name '.txt'], 'file') & exist([Y_r_smooth_name '.txt'], 'file'))
                        [r, s] = system(['/data/noel/noel6/Oct-2010/CIVET-1.2.0/progs/surface-resample2 -clobber ' surf_l_xfm ' ' Y_l_smooth_name '.txt ' Y_l_smooth_name '_rsl.txt'])
                        [r, s] = system(['/data/noel/noel6/Oct-2010/CIVET-1.2.0/progs/surface-resample2 -clobber ' surf_r_xfm ' ' Y_r_smooth_name '.txt ' Y_r_smooth_name '_rsl.txt'])
                    end
                end
            end
            
            for i = 1 : NumSubSurf
                basename = [ Prefix_cont '_' case_num_cont{k} '_white_surface_' num2str(i) ];
                
                for f = 1 : size(Feature, 2)
                    Y_l_smooth_name = [ OUTPATH case_num_cont{k} '/measurement/' basename '_left_'  num2str(NumMesh) postfix_surf '_' Modality{m} '_' Feature{f} ];
                    Y_r_smooth_name = [ OUTPATH case_num_cont{k} '/measurement/' basename '_right_'  num2str(NumMesh) postfix_surf '_' Modality{m} '_' Feature{f} ];
                    if(exist([Y_l_smooth_name '.txt'], 'file') & exist([Y_r_smooth_name '.txt'], 'file'))
                        [r, s] = system(['/data/noel/noel6/Oct-2010/CIVET-1.2.0/progs/surface-resample2 -clobber ' surf_l_xfm ' ' Y_l_smooth_name '.txt ' Y_l_smooth_name '_rsl.txt'])
                        [r, s] = system(['/data/noel/noel6/Oct-2010/CIVET-1.2.0/progs/surface-resample2 -clobber ' surf_r_xfm ' ' Y_r_smooth_name '.txt ' Y_r_smooth_name '_rsl.txt'])
                    end
                end
            end
        end
        toc
    end
    
    % patients
    for k = 1 : size(case_num_pat, 1)
        fprintf('Patient case: %s start! ', case_num_pat{k});
        tic
        
        for m = 1 : size(Modality, 2)
            Modality_temp = Modality{m};
            if(strcmp(Modality_temp, 't1_flair_ratio'))
                Modality_temp = 't1';
            end
            
            surf_l_xfm = [ OUTPATH case_num_pat{k} '/xfm/' Prefix_pat '_' case_num_pat{k} '_left_surfmap.sm' ];
            surf_r_xfm = [ OUTPATH case_num_pat{k} '/xfm/' Prefix_pat '_' case_num_pat{k} '_right_surfmap.sm' ];
            
            for i = 1 : NumIntSurf + 2
                if(i == 1)
                    basename  = [ Prefix_pat '_' case_num_pat{k} '_gray_surface' ];
                elseif(i == NumIntSurf + 2)
                    basename  = [ Prefix_pat '_' case_num_pat{k} '_white_surface' ];
                else
                    basename  = [ Prefix_pat '_' case_num_pat{k} '_intracortical_surface_' num2str(i-1) ];
                end
                
                for f = 1 : size(Feature, 2)
                    Y_l_smooth_name = [ OUTPATH case_num_pat{k} '/measurement/' basename '_left_'  num2str(NumMesh) postfix_surf '_' Modality{m} '_' Feature{f} ];
                    Y_r_smooth_name = [ OUTPATH case_num_pat{k} '/measurement/' basename '_right_'  num2str(NumMesh) postfix_surf '_' Modality{m} '_' Feature{f} ];
                    if(exist([Y_l_smooth_name '.txt'], 'file') & exist([Y_r_smooth_name '.txt'], 'file'))
                        [r, s] = system(['/data/noel/noel6/Oct-2010/CIVET-1.2.0/progs/surface-resample2 -clobber ' surf_l_xfm ' ' Y_l_smooth_name '.txt ' Y_l_smooth_name '_rsl.txt'])
                        [r, s] = system(['/data/noel/noel6/Oct-2010/CIVET-1.2.0/progs/surface-resample2 -clobber ' surf_r_xfm ' ' Y_r_smooth_name '.txt ' Y_r_smooth_name '_rsl.txt'])
                    end
                end
            end
            
            for i = 1 : NumSubSurf
                basename = [ Prefix_pat '_' case_num_pat{k} '_white_surface_' num2str(i) ];
                
                for f = 1 : size(Feature, 2)
                    Y_l_smooth_name = [ OUTPATH case_num_pat{k} '/measurement/' basename '_left_'  num2str(NumMesh) postfix_surf '_' Modality{m} '_' Feature{f} ];
                    Y_r_smooth_name = [ OUTPATH case_num_pat{k} '/measurement/' basename '_right_'  num2str(NumMesh) postfix_surf '_' Modality{m} '_' Feature{f} ];
                    if(exist([Y_l_smooth_name '.txt'], 'file') & exist([Y_r_smooth_name '.txt'], 'file'))
                        [r, s] = system(['/data/noel/noel6/Oct-2010/CIVET-1.2.0/progs/surface-resample2 -clobber ' surf_l_xfm ' ' Y_l_smooth_name '.txt ' Y_l_smooth_name '_rsl.txt'])
                        [r, s] = system(['/data/noel/noel6/Oct-2010/CIVET-1.2.0/progs/surface-resample2 -clobber ' surf_r_xfm ' ' Y_r_smooth_name '.txt ' Y_r_smooth_name '_rsl.txt'])
                    end
                end
            end
        end
        toc
    end
    
end
