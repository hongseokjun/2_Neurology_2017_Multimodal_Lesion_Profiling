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
    Kernel                    = 2
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
    
    %% feature-specific parameter tuning (weight degree for extralesion and lesion)
    aniso_smooth_parameters = { [ 1.5 3 ],          %% RI
                                [ 1.5 3 ],          %% RI_corrected
                                [ 1   2 ],          %% pg
                                [ 1   2 ],          %% pg_GM_WM
                                [ 2   3 ],          %% tg
                                [ 7   4 ],          %% ct
                                [ 7   5 ],          %% sd2
                                [ 12  12 ],         %% mc
                                [ 25  25 ],         %% FA
                                [ 10000 10000 ] };  %% MD

end

%% Smooth the data outside the given mask to be unaffected by the abnormalities in the lesion
for anisotropic_smooth = 1
    
    % Control
    % In controls, as there is no lesion mask, we just smoothed whole
    % vertices using anisotropic kernel which gives almost same result of
    % heat kernel smoothing at this time.
    for k = 1 : size(case_num_cont, 1)
        fprintf('Control case: %s start! ', case_num_cont{k});
        
        tic
        [r, s] = system(['gunzip -vf ' OUTPATH case_num_cont{k} '/surfaces/*.obj.gz']);        
        
        for m = 1 : size(Modality, 2)
            Modality_temp = Modality{m};
            if(strcmp(Modality_temp, 't1_flair_ratio'))
                Modality_temp = 't1';
            end
            
            surf_l_name = [ OUTPATH case_num_cont{k} '/surfaces/' Prefix_cont '_' case_num_cont{k} '_intracortical_surface_2_left_'  num2str(NumMesh) '_t1.obj' ];
            surf_r_name = [ OUTPATH case_num_cont{k} '/surfaces/' Prefix_cont '_' case_num_cont{k} '_intracortical_surface_2_right_' num2str(NumMesh) '_t1.obj' ];
            surf_final  = SurfStatReadSurf({ surf_l_name, surf_r_name });
            MidMask     = SurfStatMaskCut(surf_final);
            for i = 1 : NumIntSurf + 2
                if(i == 1)
                    basename  = [ Prefix_cont '_' case_num_cont{k} '_gray_surface' ];
                elseif(i == NumIntSurf + 2)
                    basename  = [ Prefix_cont '_' case_num_cont{k} '_white_surface' ];
                else
                    basename  = [ Prefix_cont '_' case_num_cont{k} '_intracortical_surface_' num2str(i-1) ];
                end
                
                surf_l_name = [ OUTPATH case_num_cont{k} '/surfaces/' basename '_left_'  num2str(NumMesh) postfix_surf '_' Modality_temp '.obj' ];
                surf_r_name = [ OUTPATH case_num_cont{k} '/surfaces/' basename '_right_' num2str(NumMesh) postfix_surf '_' Modality_temp '.obj' ];
                surf_l = SurfStatReadSurf(surf_l_name);
                surf_r = SurfStatReadSurf(surf_r_name);
                surf   = SurfStatReadSurf({ surf_l_name, surf_r_name });
                
                for f = 1 : size(Feature, 2)
                    Y_l_name = [ OUTPATH case_num_cont{k} '/measurement/' basename '_left_'  num2str(NumMesh) postfix_surf '_' Modality{m} '_' Feature{f} '.txt' ];
                    Y_r_name = [ OUTPATH case_num_cont{k} '/measurement/' basename '_right_' num2str(NumMesh) postfix_surf '_' Modality{m} '_' Feature{f} '.txt' ];
                    if(exist(Y_l_name, 'file') & exist(Y_r_name, 'file'))
                        weight_degree = aniso_smooth_parameters{f}(1);
                        Y_l = SurfStatReadData(Y_l_name);
                        Y_r = SurfStatReadData(Y_r_name);
                        Y_l(MidMask(1:40962)==0) = 0; Y_r(MidMask(40963:81924)==0) = 0;
                        Y_l_smooth = SurfStatSmoothOutsideMask( Y_l, surf_l, Kernel, weight_degree );
                        Y_r_smooth = SurfStatSmoothOutsideMask( Y_r, surf_r, Kernel, weight_degree );
                        Y_l_smooth_name = [ OUTPATH case_num_cont{k} '/measurement/' basename '_left_'  num2str(NumMesh) postfix_surf '_' Modality{m} '_' Feature{f} '_aniso_sm_' num2str(Kernel) '.txt' ];
                        Y_r_smooth_name = [ OUTPATH case_num_cont{k} '/measurement/' basename '_right_'  num2str(NumMesh) postfix_surf '_' Modality{m} '_' Feature{f} '_aniso_sm_' num2str(Kernel) '.txt' ];
                        SurfStatWriteData(Y_l_smooth_name, Y_l_smooth);
                        SurfStatWriteData(Y_r_smooth_name, Y_r_smooth);
                    end
                end
            end
            
            for i = 1 : NumSubSurf
                basename = [ Prefix_cont '_' case_num_cont{k} '_white_surface_' num2str(i) ];
                surf_l_name = [ OUTPATH case_num_cont{k} '/surfaces/' basename '_left_'  num2str(NumMesh) postfix_surf '_' Modality_temp '.obj' ];
                surf_r_name = [ OUTPATH case_num_cont{k} '/surfaces/' basename '_right_' num2str(NumMesh) postfix_surf '_' Modality_temp '.obj' ];
                if(~exist(surf_l_name, 'file') | ~exist(surf_r_name, 'file'))
                    continue;
                end
                
                surf_l = SurfStatReadSurf(surf_l_name);
                surf_r = SurfStatReadSurf(surf_r_name);
                surf   = SurfStatReadSurf({ surf_l_name, surf_r_name });
                
                for f = 1 : size(Feature, 2)
                    Y_l_name = [ OUTPATH case_num_cont{k} '/measurement/' basename '_left_'  num2str(NumMesh) postfix_surf '_' Modality{m} '_' Feature{f} '.txt' ];
                    Y_r_name = [ OUTPATH case_num_cont{k} '/measurement/' basename '_right_' num2str(NumMesh) postfix_surf '_' Modality{m} '_' Feature{f} '.txt' ];
                    if(exist(Y_l_name, 'file') & exist(Y_r_name, 'file'))
                        weight_degree = aniso_smooth_parameters{f}(1);
                        Y_l = SurfStatReadData(Y_l_name);
                        Y_r = SurfStatReadData(Y_r_name);
                        Y_l(MidMask(1:40962)==0) = 0; Y_r(MidMask(40963:81924)==0) = 0;
                        Y_l_smooth = SurfStatSmoothOutsideMask( Y_l, surf_l, Kernel, weight_degree );
                        Y_r_smooth = SurfStatSmoothOutsideMask( Y_r, surf_r, Kernel, weight_degree );
                        Y_l_smooth_name = [ OUTPATH case_num_cont{k} '/measurement/' basename '_left_'  num2str(NumMesh) postfix_surf '_' Modality{m} '_' Feature{f} '_aniso_sm_' num2str(Kernel) '.txt' ];
                        Y_r_smooth_name = [ OUTPATH case_num_cont{k} '/measurement/' basename '_right_'  num2str(NumMesh) postfix_surf '_' Modality{m} '_' Feature{f} '_aniso_sm_' num2str(Kernel) '.txt' ];
                        SurfStatWriteData(Y_l_smooth_name, Y_l_smooth);
                        SurfStatWriteData(Y_r_smooth_name, Y_r_smooth);
                    end
                end
            end
        end
        
        toc
    end
    
    % Patients
    % The areas outside and inside the mask are separately smoothed and
    % merged together at the end.
    for k = 1 : size(case_num_pat, 1)
        
        fprintf('Patient case: %s start! ', case_num_pat{k});
        
        tic
        
        [r, s] = system(['gunzip -vf ' OUTPATH case_num_pat{k} '/surfaces/*.obj.gz']);
        LesionMask_l  = SurfStatReadData([ LESIONPATH '/' Prefix_pat '_' case_num_pat{k} '_label_union_left.txt' ]);
        LesionMask_r  = SurfStatReadData([ LESIONPATH '/' Prefix_pat '_' case_num_pat{k} '_label_union_right.txt' ]);
        
        surf_l_name = [ OUTPATH case_num_pat{k} '/surfaces/' Prefix_pat '_' case_num_pat{k} '_intracortical_surface_2_left_'  num2str(NumMesh) '_t1.obj' ];
        surf_r_name = [ OUTPATH case_num_pat{k} '/surfaces/' Prefix_pat '_' case_num_pat{k} '_intracortical_surface_2_right_' num2str(NumMesh) '_t1.obj' ];
        surf_final  = SurfStatReadSurf({ surf_l_name, surf_r_name });
        MidMask     = SurfStatMaskCut(surf_final);
        
        %% Lesion area calculation ...
        temp_l = ismember(surf_final.tri, find(LesionMask_l==1));
        temp_r = ismember(surf_final.tri, find(LesionMask_r==1));
        if(sum(sum(temp_l)) > sum(sum(temp_r)))
            temp = temp_l;
        else
            temp = temp_r;
        end
        nbr_idx = find(sum(temp, 2)==3);
        lesion_area = 0;
        
        for i = 1 : size(nbr_idx, 1)
            tri_p = surf_final.tri(nbr_idx(i), :);
            lesion_area = lesion_area + tri_area(surf_final.coord(:, tri_p(1)), surf_final.coord(:, tri_p(2)), surf_final.coord(:, tri_p(3)));
        end
        
        lesion_area_set = [lesion_area_set lesion_area];
        
        kernel = Kernel;
        if(lesion_area < 100)
            %         kernel_lesion = sqrt(lesion_area)/2;
            kernel_lesion = Kernel;
        else
            kernel_lesion = Kernel;
        end
        
        for m = 1 : size(Modality, 2)
            Modality_temp = Modality{m};
            if(strcmp(Modality_temp, 't1_flair_ratio'))
                Modality_temp = 't1';
            end
            
            surf_l_name = [ OUTPATH case_num_pat{k} '/surfaces/' Prefix_pat '_' case_num_pat{k} '_intracortical_surface_2_left_'  num2str(NumMesh) '_t1.obj' ];
            surf_r_name = [ OUTPATH case_num_pat{k} '/surfaces/' Prefix_pat '_' case_num_pat{k} '_intracortical_surface_2_right_' num2str(NumMesh) '_t1.obj' ];
            surf_final  = SurfStatReadSurf({ surf_l_name, surf_r_name });
            for i = 1 : NumIntSurf + 2
                if(i == 1)
                    basename  = [ Prefix_pat '_' case_num_pat{k} '_gray_surface' ];
                elseif(i == NumIntSurf + 2)
                    basename  = [ Prefix_pat '_' case_num_pat{k} '_white_surface' ];
                else
                    basename  = [ Prefix_pat '_' case_num_pat{k} '_intracortical_surface_' num2str(i-1) ];
                end
                
                surf_l_name = [ OUTPATH case_num_pat{k} '/surfaces/' basename '_left_'  num2str(NumMesh) postfix_surf '_' Modality_temp '.obj' ];
                surf_r_name = [ OUTPATH case_num_pat{k} '/surfaces/' basename '_right_' num2str(NumMesh) postfix_surf '_' Modality_temp '.obj' ];
                if(~exist(surf_l_name, 'file') | ~exist(surf_r_name, 'file'))
                    continue;
                end
                
                surf_l = SurfStatReadSurf(surf_l_name);
                surf_r = SurfStatReadSurf(surf_r_name);
                surf   = SurfStatReadSurf({ surf_l_name, surf_r_name });
                
                for f = 1 : size(Feature, 2)
                    Y_l_name = [ OUTPATH case_num_pat{k} '/measurement/' basename '_left_'  num2str(NumMesh) postfix_surf '_' Modality{m} '_' Feature{f} '.txt' ];
                    Y_r_name = [ OUTPATH case_num_pat{k} '/measurement/' basename '_right_' num2str(NumMesh) postfix_surf '_' Modality{m} '_' Feature{f} '.txt' ];
                    if(exist(Y_l_name, 'file') & exist(Y_r_name, 'file'))
                        weight_degree = aniso_smooth_parameters{f}(1);
                        weight_degree_lesion = aniso_smooth_parameters{f}(2);
                        Y_l = SurfStatReadData(Y_l_name);
                        Y_r = SurfStatReadData(Y_r_name);
                        Y_l_org = Y_l; Y_r_org = Y_r;
                        
                        %% anisotropic smoothing outside the lesion ...
                        Y_l(LesionMask_l==1) = 0; Y_r(LesionMask_r==1) = 0;
                        Y_l(MidMask(1:40962)==0) = 0; Y_r(MidMask(40963:81924)==0) = 0;
                        Y_l_smooth_extralesion = SurfStatSmoothOutsideMask( Y_l, surf_l, Kernel, weight_degree );
                        Y_r_smooth_extralesion = SurfStatSmoothOutsideMask( Y_r, surf_r, Kernel, weight_degree );
                        
                        %% anisotropic smoothing inside the lesion ...
                        Y_l = Y_l_org; Y_r = Y_r_org;
                        Y_l(LesionMask_l==0) = 0; Y_r(LesionMask_r==0) = 0;
                        Y_l_smooth_lesion = SurfStatSmoothOutsideMask( Y_l, surf_l, kernel_lesion, weight_degree_lesion );
                        Y_r_smooth_lesion = SurfStatSmoothOutsideMask( Y_r, surf_r, kernel_lesion, weight_degree_lesion );
                        
                        %% merge them together ...
                        Y_l_smooth = Y_l_smooth_extralesion;
                        Y_r_smooth = Y_r_smooth_extralesion;
                        Y_l_smooth(LesionMask_l==1) = Y_l_smooth_lesion(LesionMask_l==1);
                        Y_r_smooth(LesionMask_r==1) = Y_r_smooth_lesion(LesionMask_r==1);
                        
                        Y_l_smooth_name = [ OUTPATH case_num_pat{k} '/measurement/' basename '_left_'  num2str(NumMesh) postfix_surf '_' Modality{m} '_' Feature{f} '_aniso_sm_' num2str(Kernel) '.txt' ];
                        Y_r_smooth_name = [ OUTPATH case_num_pat{k} '/measurement/' basename '_right_'  num2str(NumMesh) postfix_surf '_' Modality{m} '_' Feature{f} '_aniso_sm_' num2str(Kernel) '.txt' ];
                        SurfStatWriteData(Y_l_smooth_name, Y_l_smooth);
                        SurfStatWriteData(Y_r_smooth_name, Y_r_smooth);
                    end
                end
            end
            
            for i = 1 : NumSubSurf
                basename = [ Prefix_pat '_' case_num_pat{k} '_white_surface_' num2str(i) ];
                surf_l_name = [ OUTPATH case_num_pat{k} '/surfaces/' basename '_left_'  num2str(NumMesh) postfix_surf '_' Modality_temp '.obj' ];
                surf_r_name = [ OUTPATH case_num_pat{k} '/surfaces/' basename '_right_' num2str(NumMesh) postfix_surf '_' Modality_temp '.obj' ];
                if(~exist(surf_l_name, 'file') | ~exist(surf_r_name, 'file'))
                    continue;
                end
                
                surf_l = SurfStatReadSurf(surf_l_name);
                surf_r = SurfStatReadSurf(surf_r_name);
                surf   = SurfStatReadSurf({ surf_l_name, surf_r_name });
                
                for f = 1 : size(Feature, 2)
                    Y_l_name = [ OUTPATH case_num_pat{k} '/measurement/' basename '_left_'  num2str(NumMesh) postfix_surf '_' Modality{m} '_' Feature{f} '.txt' ];
                    Y_r_name = [ OUTPATH case_num_pat{k} '/measurement/' basename '_right_' num2str(NumMesh) postfix_surf '_' Modality{m} '_' Feature{f} '.txt' ];
                    if(exist(Y_l_name, 'file') & exist(Y_r_name, 'file'))
                        weight_degree = aniso_smooth_parameters{f}(1);
                        weight_degree_lesion = aniso_smooth_parameters{f}(2);
                        Y_l = SurfStatReadData(Y_l_name);
                        Y_r = SurfStatReadData(Y_r_name);
                        Y_l_org = Y_l; Y_r_org = Y_r;
                        
                        %% anisotropic smoothing outside the lesion ...
                        Y_l(LesionMask_l==1) = 0; Y_r(LesionMask_r==1) = 0;
                        Y_l(MidMask(1:40962)==0) = 0; Y_r(MidMask(40963:81924)==0) = 0;
                        Y_l_smooth_extralesion = SurfStatSmoothOutsideMask( Y_l, surf_l, Kernel, weight_degree );
                        Y_r_smooth_extralesion = SurfStatSmoothOutsideMask( Y_r, surf_r, Kernel, weight_degree );
                        
                        %% anisotropic smoothing inside the lesion ...
                        Y_l = Y_l_org; Y_r = Y_r_org;
                        Y_l(LesionMask_l==0) = 0; Y_r(LesionMask_r==0) = 0;
                        Y_l_smooth_lesion = SurfStatSmoothOutsideMask( Y_l, surf_l, kernel_lesion, weight_degree_lesion );
                        Y_r_smooth_lesion = SurfStatSmoothOutsideMask( Y_r, surf_r, kernel_lesion, weight_degree_lesion );
                        
                        %% merge them together ...
                        Y_l_smooth = Y_l_smooth_extralesion;
                        Y_r_smooth = Y_r_smooth_extralesion;
                        Y_l_smooth(LesionMask_l==1) = Y_l_smooth_lesion(LesionMask_l==1);
                        Y_r_smooth(LesionMask_r==1) = Y_r_smooth_lesion(LesionMask_r==1);
                        
                        Y_l_smooth_name = [ OUTPATH case_num_pat{k} '/measurement/' basename '_left_'  num2str(NumMesh) postfix_surf '_' Modality{m} '_' Feature{f} '_aniso_sm_' num2str(Kernel) '.txt' ];
                        Y_r_smooth_name = [ OUTPATH case_num_pat{k} '/measurement/' basename '_right_'  num2str(NumMesh) postfix_surf '_' Modality{m} '_' Feature{f} '_aniso_sm_' num2str(Kernel) '.txt' ];
                        SurfStatWriteData(Y_l_smooth_name, Y_l_smooth);
                        SurfStatWriteData(Y_r_smooth_name, Y_r_smooth);
                    end
                end
            end
        end
        
        toc
    end
    
end

%% Resample the features
for resample_features = 1
    
    %% control resampling ...
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
                    Y_l_smooth_name = [ OUTPATH case_num_cont{k} '/measurement/' basename '_left_'  num2str(NumMesh) postfix_surf '_' Modality{m} '_' Feature{f} '_aniso_sm_' num2str(Kernel) ];
                    Y_r_smooth_name = [ OUTPATH case_num_cont{k} '/measurement/' basename '_right_'  num2str(NumMesh) postfix_surf '_' Modality{m} '_' Feature{f} '_aniso_sm_' num2str(Kernel) ];
                    if(exist([Y_l_smooth_name '.txt'], 'file') & exist([Y_r_smooth_name '.txt'], 'file'))
                        [r, s] = system(['/data/noel/noel6/Oct-2010/CIVET-1.2.0/progs/surface-resample2 -clobber ' surf_l_xfm ' ' Y_l_smooth_name '.txt ' Y_l_smooth_name '_rsl.txt'])
                        [r, s] = system(['/data/noel/noel6/Oct-2010/CIVET-1.2.0/progs/surface-resample2 -clobber ' surf_r_xfm ' ' Y_r_smooth_name '.txt ' Y_r_smooth_name '_rsl.txt'])
                    end
                end
            end
            
            for i = 1 : NumSubSurf
                basename = [ Prefix_cont '_' case_num_cont{k} '_white_surface_' num2str(i) ];
                
                for f = 1 : size(Feature, 2)
                    Y_l_smooth_name = [ OUTPATH case_num_cont{k} '/measurement/' basename '_left_'  num2str(NumMesh) postfix_surf '_' Modality{m} '_' Feature{f} '_aniso_sm_' num2str(Kernel) ];
                    Y_r_smooth_name = [ OUTPATH case_num_cont{k} '/measurement/' basename '_right_'  num2str(NumMesh) postfix_surf '_' Modality{m} '_' Feature{f} '_aniso_sm_' num2str(Kernel) ];
                    if(exist([Y_l_smooth_name '.txt'], 'file') & exist([Y_r_smooth_name '.txt'], 'file'))
                        [r, s] = system(['/data/noel/noel6/Oct-2010/CIVET-1.2.0/progs/surface-resample2 -clobber ' surf_l_xfm ' ' Y_l_smooth_name '.txt ' Y_l_smooth_name '_rsl.txt'])
                        [r, s] = system(['/data/noel/noel6/Oct-2010/CIVET-1.2.0/progs/surface-resample2 -clobber ' surf_r_xfm ' ' Y_r_smooth_name '.txt ' Y_r_smooth_name '_rsl.txt'])
                    end
                end
            end
        end
        toc
    end
    
    %% patient resampling ...
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
                    Y_l_smooth_name = [ OUTPATH case_num_pat{k} '/measurement/' basename '_left_'  num2str(NumMesh) postfix_surf '_' Modality{m} '_' Feature{f} '_aniso_sm_' num2str(Kernel) ];
                    Y_r_smooth_name = [ OUTPATH case_num_pat{k} '/measurement/' basename '_right_'  num2str(NumMesh) postfix_surf '_' Modality{m} '_' Feature{f} '_aniso_sm_' num2str(Kernel) ];
                    if(exist([Y_l_smooth_name '.txt'], 'file') & exist([Y_r_smooth_name '.txt'], 'file'))
                        [r, s] = system(['/data/noel/noel6/Oct-2010/CIVET-1.2.0/progs/surface-resample2 -clobber ' surf_l_xfm ' ' Y_l_smooth_name '.txt ' Y_l_smooth_name '_rsl.txt'])
                        [r, s] = system(['/data/noel/noel6/Oct-2010/CIVET-1.2.0/progs/surface-resample2 -clobber ' surf_r_xfm ' ' Y_r_smooth_name '.txt ' Y_r_smooth_name '_rsl.txt'])
                    end
                end
            end
            
            for i = 1 : NumSubSurf
                basename = [ Prefix_pat '_' case_num_pat{k} '_white_surface_' num2str(i) ];
                
                for f = 1 : size(Feature, 2)
                    Y_l_smooth_name = [ OUTPATH case_num_pat{k} '/measurement/' basename '_left_'  num2str(NumMesh) postfix_surf '_' Modality{m} '_' Feature{f} '_aniso_sm_' num2str(Kernel) ];
                    Y_r_smooth_name = [ OUTPATH case_num_pat{k} '/measurement/' basename '_right_'  num2str(NumMesh) postfix_surf '_' Modality{m} '_' Feature{f} '_aniso_sm_' num2str(Kernel) ];
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
