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
    Cases_disease_cont        = '/local_raid/seokjun/01_project/04_IntracorticalAnalysis/01_analysis/Demographic_data_TLE_final.txt'
    Group_disease_cont        = 'TLE'
    Prefix_disease_cont       = 'TLE'
    
    Left_Right                = 'both'
    NumIntSurf                = 3
    NumSubSurf                = 5
    NumMesh                   = 81920
    SamplingSpace             = 'native'
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
    
    outlier_cases = { '306_1', '313_1', '322_1' };
    case_cont_real = 1:size(case_num_cont, 1);
    outlier_vector = zeros(size(case_num_cont, 1), 1);
    for i = 1 : size(outlier_cases, 2)
        outlier_vector = outlier_vector + strcmp(case_num_cont, outlier_cases{i});
    end
    case_num_cont(logical(outlier_vector)) = [];
    
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
    
    fid = fopen(Cases_disease_cont);
    demo = textscan(fid, '%s%f%s%s', 'Delimiter', ',', 'CollectOutput', 1);
    case_num_disease_cont       = demo{1};
    age_disease_cont            = demo{2};
    gender_disease_cont         = demo{3}(:, 1);
    seizure_lateralization  = demo{3}(:, 2);
    fclose(fid);
    
    mean_age_cont = mean(age_cont);
    std_age_cont = std(age_cont, 0);
    
    mean_age_pat = mean(age_pat);
    std_age_pat = std(age_pat, 0);
    
    [h,p,ci,stats] = ttest2(age_cont,age_pat);
    
end

%% Set up parameters related to the analysis
for setup_params = 1
    
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
        
    fMRI_surfpath    = '/local_raid/seokjun/01_project/05_NoelRest/03_Result/01_DPARSF_processing/DPARSFA_scrub/02_without_mean_sig_regout/bbr_rest/';
    lesion_label_dir = '/local_raid/seokjun/CTFCD-1.2.0_64/Lesion/lesion_surf/';
    postfix_surf     = '_native';
    Modality         = { 't1', 'flair', 'dti', 'fMRI' };
    
    geoDist_t1       = cell(size(case_num_pat, 1), 1);
    geoDist_flair    = cell(size(case_num_pat, 1), 1);
    geoDist_dti      = cell(size(case_num_pat, 1), 1);
    geoDist_fMRI     = cell(size(case_num_pat, 1), 1);
    
    geoDist_t1_cont       = cell(size(case_num_cont, 1), 1);
    geoDist_flair_cont    = cell(size(case_num_cont, 1), 1);
    geoDist_dti_cont      = cell(size(case_num_cont, 1), 1);
    geoDist_fMRI_cont     = cell(size(case_num_cont, 1), 1);
    
    geoDist_t1_disease_cont       = cell(size(case_num_disease_cont, 1), 1);
    geoDist_flair_disease_cont    = cell(size(case_num_disease_cont, 1), 1);
    geoDist_dti_disease_cont      = cell(size(case_num_disease_cont, 1), 1);
    geoDist_fMRI_disease_cont     = cell(size(case_num_disease_cont, 1), 1);
        
    surf_ind         = cell(size(case_num_pat, 1), 1);
    surf_ind_t1      = cell(size(case_num_pat, 1), 1);
    surf_ind_flair   = cell(size(case_num_pat, 1), 1);
    surf_ind_dti     = cell(size(case_num_pat, 1), 1);
    surf_ind_fMRI    = cell(size(case_num_pat, 1), 1);
    
    surf_ind_cont         = cell(size(case_num_cont, 1), 1);
    surf_ind_t1_cont      = cell(size(case_num_cont, 1), 1);
    surf_ind_flair_cont   = cell(size(case_num_cont, 1), 1);
    surf_ind_dti_cont     = cell(size(case_num_cont, 1), 1);
    surf_ind_fMRI_cont    = cell(size(case_num_cont, 1), 1);
        
    surf_ind_t1_disease_cont      = cell(size(case_num_disease_cont, 1), 1);
    surf_ind_flair_disease_cont   = cell(size(case_num_disease_cont, 1), 1);
    surf_ind_dti_disease_cont     = cell(size(case_num_disease_cont, 1), 1);
    surf_ind_fMRI_disease_cont    = cell(size(case_num_disease_cont, 1), 1);
    
    surf   = SurfStatReadSurf({ '/local_raid/seokjun/01_project/04_IntracorticalAnalysis/03_result/002/surfaces/mcd_002_gray_surface_left_81920_native_t1.obj', ...
                                '/local_raid/seokjun/01_project/04_IntracorticalAnalysis/03_result/002/surfaces/mcd_002_gray_surface_right_81920_native_t1.obj' });
    surf = surfGetNeighborsHong(surf);
    nbr = surf.nbr';

end

%% Compute geodesic distance for the lesion
for geodesic_dist_lesion = 1
    
    for patient_comp = 1
        
        OUTPATH          = '/local_raid/seokjun/01_project/04_IntracorticalAnalysis/03_result/';        
        case_num_pat_org = case_num_pat;
        
        % 1) The original method
        for j = 1 : size(case_num_pat, 1)
            idx = find(strcmp(case_num_pat, case_num_pat{j}));
            if(exist([ lesion_label_dir '/' Prefix_pat '_' case_num_pat{j} '_label_union_left_rsl.txt' ], 'file'))
                lesion_label_data = SurfStatReadData( { [ lesion_label_dir '/' Prefix_pat '_' case_num_pat{j} '_label_union_left_rsl.txt' ], ...
                    [ lesion_label_dir '/' Prefix_pat '_' case_num_pat{j} '_label_union_right_rsl.txt' ] } );
                for m = 1 : size(Modality, 2)
                    
                    for T1_flair_dti = 1
                        if(strcmp(Modality{m}, 't1') | strcmp(Modality{m}, 'flair') | strcmp(Modality{m}, 'dti'))
                            if(strcmp(case_num_pat{j}, '080_1') & strcmp(Modality{m}, 'dti'))
                                case_num_pat{j} = '080_2';
                            end
                            eval([ 'geoDist_' Modality{m} '{' num2str(j) '}.surf = zeros(NumIntSurf + 2 + NumSubSurf, 81924);' ]);
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
                                
                                surf_l_name = [ OUTPATH case_num_pat{j} '/surfaces/' basename '_left_'  num2str(NumMesh) postfix_surf '_' Modality{m} '.obj' ];
                                surf_r_name = [ OUTPATH case_num_pat{j} '/surfaces/' basename '_right_' num2str(NumMesh) postfix_surf '_' Modality{m} '.obj' ];
                                surf_l_rsl_name = [ OUTPATH case_num_pat{j} '/surfaces/' basename '_left_'  num2str(NumMesh) postfix_surf '_' Modality{m} '_rsl.obj' ];
                                surf_r_rsl_name = [ OUTPATH case_num_pat{j} '/surfaces/' basename '_right_' num2str(NumMesh) postfix_surf '_' Modality{m} '_rsl.obj' ];
                                xfm_l_name  = [ OUTPATH case_num_pat{j} '/xfm/' Prefix_pat '_' case_num_pat{j} '_left_surfmap.sm' ];
                                xfm_r_name  = [ OUTPATH case_num_pat{j} '/xfm/' Prefix_pat '_' case_num_pat{j} '_right_surfmap.sm' ];
                                
                                if(exist(surf_l_rsl_name, 'file') & exist(surf_r_rsl_name, 'file'))
                                    disp([ basename '_' Modality{m} ' exists.' ]);
                                    surf   = SurfStatReadSurf({ surf_l_rsl_name, surf_r_rsl_name });
                                    surf_ind{j}.surf{i}.tri = surf.tri;
                                    surf_ind{j}.surf{i}.coord = surf.coord;
                                    surf.nbr = nbr;
                                    eval([ 'surf_ind_' Modality{m} '{' num2str(j) '}.surf{' num2str(i) '} = surf_ind{' num2str(j) '}.surf{' num2str(i) '};' ]);
                                    eval([ 'geoDist_' Modality{m} '{' num2str(j) '}.surf(' num2str(i) ', :) = surfGeoDistHong(surf, lesion_label_data); ']);
                                    continue;
                                end
                                
                                if(~exist(surf_l_name, 'file') | ~exist(surf_r_name, 'file'))
                                    continue;
                                end
                                disp([ basename '_' Modality{m} ' is resampled now!' ]);
                                [r, s] = system(['sphere_resample_obj ' surf_l_name ' ' xfm_l_name ' ' surf_l_rsl_name ]);
                                [r, s] = system(['sphere_resample_obj ' surf_r_name ' ' xfm_r_name ' ' surf_r_rsl_name ]);
                                
                                surf   = SurfStatReadSurf({ surf_l_rsl_name, surf_r_rsl_name });
                                surf_ind{j}.surf{i}.tri = surf.tri;
                                surf_ind{j}.surf{i}.coord = surf.coord;
                                surf.nbr = nbr;
                                eval([ 'surf_ind_' Modality{m} '{' num2str(j) '}.surf{' num2str(i) '} = surf_ind{' num2str(j) '}.surf{' num2str(i) '};' ]);
                                eval([ 'geoDist_' Modality{m} '{' num2str(j) '}.surf(' num2str(i) ', :) = surfGeoDistHong(surf, lesion_label_data); ']);
                            end
                            for i = 1 : NumSubSurf
                                basename = [ Prefix_pat '_' case_num_pat{j} '_white_surface_' num2str(i) ];
                                temp_str = num2str(i);
                                
                                surf_l_name = [ OUTPATH case_num_pat{j} '/surfaces/' basename '_left_'  num2str(NumMesh) postfix_surf '_' Modality{m} '.obj' ];
                                surf_r_name = [ OUTPATH case_num_pat{j} '/surfaces/' basename '_right_' num2str(NumMesh) postfix_surf '_' Modality{m} '.obj' ];
                                surf_l_rsl_name = [ OUTPATH case_num_pat{j} '/surfaces/' basename '_left_'  num2str(NumMesh) postfix_surf '_' Modality{m} '_rsl.obj' ];
                                surf_r_rsl_name = [ OUTPATH case_num_pat{j} '/surfaces/' basename '_right_' num2str(NumMesh) postfix_surf '_' Modality{m} '_rsl.obj' ];
                                xfm_l_name  = [ OUTPATH case_num_pat{j} '/xfm/' Prefix_pat '_' case_num_pat{j} '_left_surfmap.sm' ];
                                xfm_r_name  = [ OUTPATH case_num_pat{j} '/xfm/' Prefix_pat '_' case_num_pat{j} '_right_surfmap.sm' ];
                                
                                if(exist(surf_l_rsl_name, 'file') & exist(surf_r_rsl_name, 'file'))
                                    disp([ basename '_' Modality{m} ' exists.' ]);
                                    surf   = SurfStatReadSurf({ surf_l_rsl_name, surf_r_rsl_name });
                                    surf_ind{j}.surf{i+5}.tri = surf.tri;
                                    surf_ind{j}.surf{i+5}.coord = surf.coord;
                                    surf.nbr = nbr;
                                    eval([ 'surf_ind_' Modality{m} '{' num2str(j) '}.surf{' num2str(i+5) '} = surf_ind{' num2str(j) '}.surf{' num2str(i+5) '};' ]);
                                    eval([ 'geoDist_' Modality{m} '{' num2str(j) '}.surf(' num2str(i+5) ', :) = surfGeoDistHong(surf, lesion_label_data); ']);
                                    continue;
                                end
                                
                                if(~exist(surf_l_name, 'file') | ~exist(surf_r_name, 'file'))
                                    continue;
                                end
                                disp([ basename '_' Modality{m} ' is resampled now!' ]);
                                [r, s] = system(['sphere_resample_obj ' surf_l_name ' ' xfm_l_name ' ' surf_l_rsl_name ]);
                                [r, s] = system(['sphere_resample_obj ' surf_r_name ' ' xfm_r_name ' ' surf_r_rsl_name ]);
                                
                                surf   = SurfStatReadSurf({ surf_l_rsl_name, surf_r_rsl_name });
                                surf_ind{j}.surf{i+5}.tri = surf.tri;
                                surf_ind{j}.surf{i+5}.coord = surf.coord;
                                surf.nbr = nbr;
                                eval([ 'surf_ind_' Modality{m} '{' num2str(j) '}.surf{' num2str(i+5) '} = surf_ind{' num2str(j) '}.surf{' num2str(i+5) '};' ]);
                                eval([ 'geoDist_' Modality{m} '{' num2str(j) '}.surf(' num2str(i+5) ', :) = surfGeoDistHong(surf, lesion_label_data); ']);
                            end
                        end
                    end
                    
                    for fMRI = 1
                        
                        if(strcmp(Modality{m}, 'fMRI'))
                            if(strcmp(case_num_pat{j}, '080_1') & strcmp(Modality{m}, 'fMRI'))
                                case_num_pat{j} = '080_2';
                            end
                            eval([ 'geoDist_' Modality{m} '{' num2str(j) '}.surf = zeros(1, 81924);' ]);
                            surf_ind{j}.surf = cell(1, 1);
                            
                            surf_l_name = [ fMRI_surfpath '/' Prefix_pat '_' case_num_pat{j} '/' Prefix_pat '_' case_num_pat{j} '_mid_surface_left_' num2str(NumMesh)  '_rest_bbr.obj' ];
                            surf_r_name = [ fMRI_surfpath '/' Prefix_pat '_' case_num_pat{j} '/' Prefix_pat '_' case_num_pat{j} '_mid_surface_right_' num2str(NumMesh) '_rest_bbr.obj' ];
                            surf_l_rsl_name = [ fMRI_surfpath '/' Prefix_pat '_' case_num_pat{j} '/' Prefix_pat '_' case_num_pat{j} '_mid_surface_left_' num2str(NumMesh)  '_rest_bbr_rsl.obj' ];
                            surf_r_rsl_name = [ fMRI_surfpath '/' Prefix_pat '_' case_num_pat{j} '/' Prefix_pat '_' case_num_pat{j} '_mid_surface_right_' num2str(NumMesh) '_rest_bbr_rsl.obj' ];
                            
                            xfm_l_name  = [ OUTPATH case_num_pat{j} '/xfm/' Prefix_pat '_' case_num_pat{j} '_left_surfmap.sm' ];
                            xfm_r_name  = [ OUTPATH case_num_pat{j} '/xfm/' Prefix_pat '_' case_num_pat{j} '_right_surfmap.sm' ];
                            
                            if(exist(surf_l_rsl_name, 'file') & exist(surf_r_rsl_name, 'file'))
                                disp([ fMRI_surfpath '/' Prefix_pat '_' case_num_pat{j} '/' Prefix_pat '_' case_num_pat{j} '_mid_surface_rest_bbr_rsl exists.' ]);
                                surf   = SurfStatReadSurf({ surf_l_rsl_name, surf_r_rsl_name });
                                surf_ind{j}.surf{1}.tri = surf.tri;
                                surf_ind{j}.surf{1}.coord = surf.coord;
                                surf.nbr = nbr;
                                eval([ 'surf_ind_' Modality{m} '{' num2str(j) '}.surf{1} = surf_ind{' num2str(j) '}.surf{1};' ]);
                                eval([ 'geoDist_' Modality{m} '{' num2str(j) '}.surf(1, :) = surfGeoDistHong(surf, lesion_label_data); ']);
                                continue;
                            end
                            
                            if(~exist(surf_l_name, 'file') | ~exist(surf_r_name, 'file'))
                                continue;
                            end
                            disp([ fMRI_surfpath '/' Prefix_pat '_' case_num_pat{j} '/' Prefix_pat '_' case_num_pat{j} '_mid_surface_rest_bbr is resampled now!' ]);
                            [r, s] = system(['sphere_resample_obj ' surf_l_name ' ' xfm_l_name ' ' surf_l_rsl_name ]);
                            [r, s] = system(['sphere_resample_obj ' surf_r_name ' ' xfm_r_name ' ' surf_r_rsl_name ]);
                            
                            surf   = SurfStatReadSurf({ surf_l_rsl_name, surf_r_rsl_name });
                            surf_ind{j}.surf{1}.tri = surf.tri;
                            surf_ind{j}.surf{1}.coord = surf.coord;
                            surf.nbr = nbr;
                            eval([ 'surf_ind_' Modality{m} '{' num2str(j) '}.surf{1} = surf_ind{' num2str(j) '}.surf{1};' ]);
                            eval([ 'geoDist_' Modality{m} '{' num2str(j) '}.surf(1, :) = surfGeoDistHong(surf, lesion_label_data); ']);
                            
                        end
                    end
                    
                end
                disp(['case : ' case_num_pat{j} ' done']);
            end
        end
        
        case_num_pat = case_num_pat_org;
        
        % 2) Modified geodesic distance
        % Problem: An originally calculated geodesic distance on subcortical white matter surface has an intrinsic limitation
        % at the narrow white matter region where two surfaces are jammed together very closely.
        % It's because although the original geodesic distance is gradually growing along the surface as far from the lesion,
        % sampling the value is being done still in same area.
        % Solution: 1) We apply a small sphere kernel at every vertex to detect such jammed surface areas.
        %              If more than two vertices are included in the sphere, we assume the surfaces at that location are likely jammed each other.
        %           2) We assign minimal geodesic distance to all detected vertices.
        sp_kernel = 2;
        surf = surfGetNeighborsHong(surf);
        nbr = surf.nbr;
        load('individual_surface_multimodal.mat');
        geoDist_t1_new = geoDist_t1;
        geoDist_flair_new = geoDist_flair;
        geoDist_dti_new = geoDist_dti;
        geoDist_fMRI_new = geoDist_fMRI;
        for j = 1 : size(case_num_pat, 1)
            idx = find(strcmp(case_num_pat, case_num_pat{j}));
            for m = 1 : size(Modality, 2)
                eval([ 'geoDist  = geoDist_' Modality{m} '{' num2str(j) '};'  ]);
                eval([ 'surf_ind = surf_ind_' Modality{m} '{' num2str(j) '};'  ]);
                if(strcmp(case_num_pat{j}, '080_1') & (strcmp(Modality{m}, 'dti') | strcmp(Modality{m}, 'fMRI')))
                    case_num_pat{j} = '080_2';
                end
                
                if(strcmp(Modality{m}, 't1') | strcmp(Modality{m}, 'flair'))
                    NumSubSurf_temp = 3;
                elseif(strcmp(Modality{m}, 'dti'))
                    NumSubSurf_temp = 5;
                elseif(strcmp(Modality{m}, 'fMRI'))
                    NumSubSurf_temp = 0;
                end
                
                for i = 1 : NumSubSurf_temp
                    jammed_vertices = [];
                    idx_vertex = find(geoDist.surf(i + 5, :) < 50 & geoDist.surf(i + 5, :) > 0); %% limited the circal range where we detect a jammed surface area within 50mm
                    geoDist_temp = geoDist.surf(i + 5, :);
                    for iv = 1 : size(idx_vertex, 2)
                        center = surf_ind.surf{i + 5}.coord(:, idx_vertex(iv));                                                             % for those vertices
                        distance = sqrt(sum((surf_ind.surf{i + 5}.coord - kron(ones(1, size(surf_ind.surf{i + 5}.coord, 2)), center)).^2)); % calculate a eucliean distance for all the 80k vertices at once
                        jammed_vertices_temp = find(distance < sp_kernel);                                                                  % detect the only vertices that are distant within 2mm from the center
                        [a b] = intersect(jammed_vertices_temp, nbr(idx_vertex(iv), :));                                                    % if they are just first order neighbourhoods, then
                        if(~isempty(b))
                            jammed_vertices_temp(b) = [];                                                                                   % then don't consider as a candidate so exclude them in the cue
                        end
                        jammed_vertices_temp(jammed_vertices_temp == idx_vertex(iv)) = [];                                                  % also exclude the center itself
                        geoDist_temp(jammed_vertices_temp) = min(geoDist.surf(i + 5, [idx_vertex(iv) jammed_vertices_temp]));
                        jammed_vertices = [ jammed_vertices jammed_vertices_temp ];
                    end
                    eval([ 'geoDist_' Modality{m} '_new{' num2str(j) '}.surf(' num2str(i + 5) ', :) = geoDist_temp;'  ]);
                end
            end
            disp(['case : ' case_num_pat{j} ' done']);
        end
        case_num_pat = case_num_pat_org;
        
    end
    
    for control_comp = 1
        
        for healthy_control = 1
            
            OUTPATH          = '/local_raid/seokjun/01_project/04_IntracorticalAnalysis/03_result/';            
            fMRI_surfpath    = '/local_raid/seokjun/01_project/05_NoelRest/03_Result/01_DPARSF_processing/DPARSFA_scrub/02_without_mean_sig_regout/bbr_rest/';
            case_num_cont_org = case_num_cont;
            
            for j = 1 : size(case_num_pat, 1)
                disp(['patient: ', num2str(j)]);
                idx = find(strcmp(case_num_pat, case_num_pat{j}));
                if(exist([ lesion_label_dir '/' Prefix_pat '_' case_num_pat{j} '_label_union_left_rsl.txt' ], 'file'))
                    lesion_label_data_temp = SurfStatReadData( { [ lesion_label_dir '/' Prefix_pat '_' case_num_pat{j} '_label_union_left_rsl.txt' ], ...
                        [ lesion_label_dir '/' Prefix_pat '_' case_num_pat{j} '_label_union_right_rsl.txt' ] } );
                    eval([ 'lesion_label_data(', num2str(j), ', :) = lesion_label_data_temp;']);
                end
            end
            
            parpool(24);
            
            % 1) The original method
            for c = 1 : size(case_num_cont, 1)
                
                for m = 1 : size(Modality, 2)
                    
                    for T1_flair_dti = 1
                        
                        if(strcmp(Modality{m}, 't1') | strcmp(Modality{m}, 'flair'))
                            eval([ 'geoDist_' Modality{m} '_cont{' num2str(c) '}.surf = zeros(NumIntSurf + NumSubSurf, 81924,  size(case_num_pat, 1));' ]);
                        end
                        if(strcmp(Modality{m}, 'dti'))
                            eval([ 'geoDist_' Modality{m} '_cont{' num2str(c) '}.surf = zeros(NumIntSurf + 2 + NumSubSurf, 81924,  size(case_num_pat, 1));' ]);
                        end
                        
                        surf_ind_cont{c}.surf = cell(NumIntSurf + 2 + NumSubSurf, 1);
                        
                        for i = 1 : NumIntSurf + 2
                            if(i == 1)
                                basename  = [ Prefix_cont '_' case_num_cont{c} '_gray_surface' ];
                                temp_str = 'gm';
                            elseif(i == NumIntSurf + 2)
                                basename  = [ Prefix_cont '_' case_num_cont{c} '_white_surface' ];
                                temp_str = 'wm';
                            else
                                basename  = [ Prefix_cont '_' case_num_cont{c} '_intracortical_surface_' num2str(i-1) ];
                                temp_str = num2str(i-1);
                            end
                            
                            surf_l_name = [ OUTPATH case_num_cont{c} '/surfaces/' basename '_left_'  num2str(NumMesh) postfix_surf '_' Modality{m} '.obj' ];
                            surf_r_name = [ OUTPATH case_num_cont{c} '/surfaces/' basename '_right_' num2str(NumMesh) postfix_surf '_' Modality{m} '.obj' ];
                            surf_l_rsl_name = [ OUTPATH case_num_cont{c} '/surfaces/' basename '_left_'  num2str(NumMesh) postfix_surf '_' Modality{m} '_rsl.obj' ];
                            surf_r_rsl_name = [ OUTPATH case_num_cont{c} '/surfaces/' basename '_right_' num2str(NumMesh) postfix_surf '_' Modality{m} '_rsl.obj' ];
                            xfm_l_name  = [ OUTPATH case_num_cont{c} '/xfm/' Prefix_cont '_' case_num_cont{c} '_left_surfmap.sm' ];
                            xfm_r_name  = [ OUTPATH case_num_cont{c} '/xfm/' Prefix_cont '_' case_num_cont{c} '_right_surfmap.sm' ];
                            
                            if(exist(surf_l_rsl_name, 'file') & exist(surf_r_rsl_name, 'file'))
                                disp([ basename '_' Modality{m} ' exists.' ]);
                                surf   = SurfStatReadSurf({ surf_l_rsl_name, surf_r_rsl_name });
                                surf_ind_cont{c}.surf{i}.tri = surf.tri;
                                surf_ind_cont{c}.surf{i}.coord = surf.coord;
                                surf.nbr = nbr;
                                eval([ 'surf_ind_' Modality{m} '_cont{' num2str(c) '}.surf{' num2str(i) '}    = surf_ind_cont{' num2str(c) '}.surf{' num2str(i) '};' ]);
                                geoDist_cont_temp = cell(size(case_num_pat, 1), 1);
                                parfor j = 1 : size(case_num_pat, 1)
                                    disp(['control: ', num2str(c), ', modality: ', num2str(m), ' surface:, ', num2str(i), ', patient: ', num2str(j)]);
                                    geoDist_cont_temp{j} = surfGeoDistHong(surf, lesion_label_data(j, :));
                                end
                                geoDist_cont_temp = cell2mat(geoDist_cont_temp);
                                eval([ 'geoDist_'  Modality{m} '_cont{' num2str(c) '}.surf(' num2str(i) ', :, :) = geoDist_cont_temp''; ' ]);
                                continue;
                            end
                            
                            if(~exist(surf_l_name, 'file') | ~exist(surf_r_name, 'file'))
                                continue;
                            end
                            
                            disp([ basename '_' Modality{m} ' is resampled now!' ]);
                            [r, s] = system(['sphere_resample_obj ' surf_l_name ' ' xfm_l_name ' ' surf_l_rsl_name ]);
                            [r, s] = system(['sphere_resample_obj ' surf_r_name ' ' xfm_r_name ' ' surf_r_rsl_name ]);
                            
                            surf   = SurfStatReadSurf({ surf_l_rsl_name, surf_r_rsl_name });
                            surf_ind_cont{c}.surf{i}.tri   = surf.tri;
                            surf_ind_cont{c}.surf{i}.coord = surf.coord;
                            surf.nbr = nbr;
                            eval([ 'surf_ind_' Modality{m} '_cont{' num2str(c) '}.surf{' num2str(i) '}    = surf_ind_cont{' num2str(c) '}.surf{' num2str(i) '};' ]);
                            geoDist_cont_temp = cell(size(case_num_pat, 1), 1);
                            parfor j = 1 : size(case_num_pat, 1)
                                disp(['control: ', num2str(c), ', modality: ', num2str(m), ' surface:, ', num2str(i), ', patient: ', num2str(j)]);
                                geoDist_cont_temp{j} = surfGeoDistHong(surf, lesion_label_data(j, :));
                            end
                            geoDist_cont_temp = cell2mat(geoDist_cont_temp);
                            eval([ 'geoDist_'  Modality{m} '_cont{' num2str(c) '}.surf(' num2str(i) ', :, :) = geoDist_cont_temp''; ' ]);
                            
                        end
                        
                        for i = 1 : NumSubSurf
                            basename = [ Prefix_cont '_' case_num_cont{c} '_white_surface_' num2str(i) ];
                            temp_str = num2str(i);
                            
                            surf_l_name = [ OUTPATH case_num_cont{c} '/surfaces/' basename '_left_'  num2str(NumMesh) postfix_surf '_' Modality{m} '.obj' ];
                            surf_r_name = [ OUTPATH case_num_cont{c} '/surfaces/' basename '_right_' num2str(NumMesh) postfix_surf '_' Modality{m} '.obj' ];
                            surf_l_rsl_name = [ OUTPATH case_num_cont{c} '/surfaces/' basename '_left_'  num2str(NumMesh) postfix_surf '_' Modality{m} '_rsl.obj' ];
                            surf_r_rsl_name = [ OUTPATH case_num_cont{c} '/surfaces/' basename '_right_' num2str(NumMesh) postfix_surf '_' Modality{m} '_rsl.obj' ];
                            xfm_l_name  = [ OUTPATH case_num_cont{c} '/xfm/' Prefix_cont '_' case_num_cont{c} '_left_surfmap.sm' ];
                            xfm_r_name  = [ OUTPATH case_num_cont{c} '/xfm/' Prefix_cont '_' case_num_cont{c} '_right_surfmap.sm' ];
                            
                            if(exist(surf_l_rsl_name, 'file') & exist(surf_r_rsl_name, 'file'))
                                disp([ basename '_' Modality{m} ' exists.' ]);
                                surf   = SurfStatReadSurf({ surf_l_rsl_name, surf_r_rsl_name });
                                surf_ind_cont{c}.surf{i+5}.tri = surf.tri;
                                surf_ind_cont{c}.surf{i+5}.coord = surf.coord;
                                surf.nbr = nbr;
                                eval([ 'surf_ind_' Modality{m} '_cont{' num2str(c) '}.surf{' num2str(i+5) '} = surf_ind_cont{' num2str(c) '}.surf{' num2str(i+5) '};' ]);
                                geoDist_cont_temp = cell(size(case_num_pat, 1), 1);
                                parfor j = 1 : size(case_num_pat, 1)
                                    disp(['control: ', num2str(c), ', modality: ', num2str(m), ' surface:, ', num2str(i+5), ', patient: ', num2str(j)]);
                                    geoDist_cont_temp{j} = surfGeoDistHong(surf, lesion_label_data(j, :));
                                end
                                geoDist_cont_temp = cell2mat(geoDist_cont_temp);
                                eval([ 'geoDist_'  Modality{m} '_cont{' num2str(c) '}.surf(' num2str(i+5) ', :, :) = geoDist_cont_temp''; ' ]);
                                continue;
                            end
                            
                            if(~exist(surf_l_name, 'file') | ~exist(surf_r_name, 'file'))
                                continue;
                            end
                            disp([ basename '_' Modality{m} ' is resampled now!' ]);
                            [r, s] = system(['sphere_resample_obj ' surf_l_name ' ' xfm_l_name ' ' surf_l_rsl_name ]);
                            [r, s] = system(['sphere_resample_obj ' surf_r_name ' ' xfm_r_name ' ' surf_r_rsl_name ]);
                            
                            surf   = SurfStatReadSurf({ surf_l_rsl_name, surf_r_rsl_name });
                            surf_ind_cont{c}.surf{i+5}.tri = surf.tri;
                            surf_ind_cont{c}.surf{i+5}.coord = surf.coord;
                            surf.nbr = nbr;
                            eval([ 'surf_ind_' Modality{m} '_cont{' num2str(c) '}.surf{' num2str(i+5) '} = surf_ind_cont{' num2str(c) '}.surf{' num2str(i+5) '};' ]);
                            geoDist_cont_temp = cell(size(case_num_pat, 1), 1);
                            parfor j = 1 : size(case_num_pat, 1)
                                disp(['control: ', num2str(c), ', modality: ', num2str(m), ' surface:, ', num2str(i+5), ', patient: ', num2str(j)]);
                                geoDist_cont_temp{j} = surfGeoDistHong(surf, lesion_label_data(j, :));
                            end
                            geoDist_cont_temp = cell2mat(geoDist_cont_temp);
                            eval([ 'geoDist_'  Modality{m} '_cont{' num2str(c) '}.surf(' num2str(i+5) ', :, :) = geoDist_cont_temp''; ' ]);
                        end
                        
                    end
                    
                    for fMRI = 1
                        
                        if(strcmp(Modality{m}, 'fMRI'))
                            
                            eval([ 'geoDist_' Modality{m} '_cont{' num2str(c) '}.surf = zeros(size(case_num_pat, 1), 81924);' ]);
                            
                            surf_ind_cont{c}.surf = cell(1, 1);
                            
                            surf_l_name = [ fMRI_surfpath '/' Prefix_cont '_' case_num_cont{c} '/' Prefix_cont '_' case_num_cont{c} '_mid_surface_left_' num2str(NumMesh)  '_rest_bbr.obj' ];
                            surf_r_name = [ fMRI_surfpath '/' Prefix_cont '_' case_num_cont{c} '/' Prefix_cont '_' case_num_cont{c} '_mid_surface_right_' num2str(NumMesh) '_rest_bbr.obj' ];
                            surf_l_rsl_name = [ fMRI_surfpath '/' Prefix_cont '_' case_num_cont{c} '/' Prefix_cont '_' case_num_cont{c} '_mid_surface_left_' num2str(NumMesh)  '_rest_bbr_rsl.obj' ];
                            surf_r_rsl_name = [ fMRI_surfpath '/' Prefix_cont '_' case_num_cont{c} '/' Prefix_cont '_' case_num_cont{c} '_mid_surface_right_' num2str(NumMesh) '_rest_bbr_rsl.obj' ];
                            
                            xfm_l_name  = [ OUTPATH case_num_cont{c} '/xfm/' Prefix_cont '_' case_num_cont{c} '_left_surfmap.sm' ];
                            xfm_r_name  = [ OUTPATH case_num_cont{c} '/xfm/' Prefix_cont '_' case_num_cont{c} '_right_surfmap.sm' ];
                            
                            if(exist(surf_l_rsl_name, 'file') & exist(surf_r_rsl_name, 'file'))
                                disp([ Prefix_cont '_' case_num_cont{c} '_mid_surface_rest_bbr_rsl exists.' ]);
                                surf   = SurfStatReadSurf({ surf_l_rsl_name, surf_r_rsl_name });
                                surf_ind_cont{c}.surf{1}.tri = surf.tri;
                                surf_ind_cont{c}.surf{1}.coord = surf.coord;
                                surf.nbr = nbr;
                                eval([ 'surf_ind_' Modality{m} '_cont{' num2str(c) '}.surf{1} = surf_ind_cont{' num2str(c) '}.surf{1};' ]);
                                geoDist_cont_temp = cell(size(case_num_pat, 1), 1);
                                parfor j = 1 : size(case_num_pat, 1)
                                    disp(['control: ', num2str(c), ', modality: ', num2str(m), ' surface:, mid, patient: ', num2str(j)]);
                                    geoDist_cont_temp{j} = surfGeoDistHong(surf, lesion_label_data(j, :));
                                end
                                geoDist_cont_temp = cell2mat(geoDist_cont_temp);
                                eval([ 'geoDist_'  Modality{m} '_cont{' num2str(c) '}.surf = geoDist_cont_temp; ' ]);
                                continue;
                            end
                            
                            if(~exist(surf_l_name, 'file') | ~exist(surf_r_name, 'file'))
                                continue;
                            end
                            
                            disp([ Prefix_cont '_' case_num_cont{c} '_mid_surface_rest_bbr is resampled now!' ]);
                            [r, s] = system(['sphere_resample_obj ' surf_l_name ' ' xfm_l_name ' ' surf_l_rsl_name ]);
                            [r, s] = system(['sphere_resample_obj ' surf_r_name ' ' xfm_r_name ' ' surf_r_rsl_name ]);
                            
                            surf   = SurfStatReadSurf({ surf_l_rsl_name, surf_r_rsl_name });
                            surf_ind_cont{c}.surf{1}.tri = surf.tri;
                            surf_ind_cont{c}.surf{1}.coord = surf.coord;
                            surf.nbr = nbr;
                            eval([ 'surf_ind_' Modality{m} '_cont{' num2str(c) '}.surf{1} = surf_ind_cont{' num2str(c) '}.surf{1};' ]);
                            geoDist_cont_temp = cell(size(case_num_pat, 1), 1);
                            parfor j = 1 : size(case_num_pat, 1)
                                disp(['control: ', num2str(c), ', modality: ', num2str(m), ' surface:, mid, patient: ', num2str(j)]);
                                geoDist_cont_temp{j} = surfGeoDistHong(surf, lesion_label_data(j, :));
                            end
                            geoDist_cont_temp = cell2mat(geoDist_cont_temp);
                            eval([ 'geoDist_'  Modality{m} '_cont{' num2str(c) '}.surf = geoDist_cont_temp; ' ]);
                            
                        end
                        
                    end
                    
                end
                
                disp(['case : ' case_num_cont{c} ' done']);
                
            end
            
            delete(gcp);
            
            case_num_cont = case_num_cont_org;
            
            % 2) Modified geodesic distance
            % Problem: An originally calculated geodesic distance on subcortical white matter surface has an intrinsic limitation
            % at the narrow white matter region where two surfaces are jammed together very closely.
            % It's because although the original geodesic distance is gradually growing along the surface as far from the lesion,
            % sampling the value is being done still in same area.
            % Solution: 1) We apply a small sphere kernel at every vertex to detect such jammed surface areas.
            %              If more than two vertices are included in the sphere, we assume the surfaces at that location are likely jammed each other.
            %           2) We assign minimal geodesic distance to all detected vertices.
            sp_kernel = 2;
            surf = surfGetNeighborsHong(surf);
            nbr = surf.nbr;
            load('individual_surface_multimodal_cont.mat');
            
            parpool(24);
            
            geoDist_t1_cont_new = geoDist_t1_cont;
            geoDist_flair_cont_new = geoDist_flair_cont;
            geoDist_dti_cont_new = geoDist_dti_cont;
            geoDist_fMRI_cont_new = geoDist_fMRI_cont;
            for j = 1 : size(case_num_cont, 1)
                
                for m = 1 : size(Modality, 2)
                    
                    eval([ 'geoDist  = geoDist_' Modality{m} '_cont{' num2str(j) '};'  ]);
                    eval([ 'surf_ind = surf_ind_' Modality{m} '_cont{' num2str(j) '};'  ]);
                    
                    if(strcmp(Modality{m}, 't1') | strcmp(Modality{m}, 'flair'))
                        NumSubSurf_temp = 3;
                    elseif(strcmp(Modality{m}, 'dti'))
                        NumSubSurf_temp = 5;
                    elseif(strcmp(Modality{m}, 'fMRI'))
                        NumSubSurf_temp = 0;
                    end
                    
                    for i = 1 : NumSubSurf_temp
                        
                        tic
                        geoDist_temp2 = cell(size(case_num_pat, 1), 1);
                        parfor p = 1 : size(case_num_pat, 1)
                            disp(['cont : ', case_num_cont{j}, ', mod : ', Modality{m}, ', surf : ', num2str(i+5) ', pat : ', case_num_pat{p}]);
                            idx_vertex = find(geoDist.surf(i + 5, :, p) < 25 & geoDist.surf(i + 5, :, p) > 0);                                      % limited the circal range where we detect a jammed surface area within 50mm
                            jammed_vertices = [];
                            geoDist_temp = geoDist.surf(i + 5, :, p);
                            for iv = 1 : size(idx_vertex, 2)
                                center = surf_ind.surf{i + 5}.coord(:, idx_vertex(iv));                                                             % for those vertices
                                distance = sqrt(sum((surf_ind.surf{i + 5}.coord - kron(ones(1, size(surf_ind.surf{i + 5}.coord, 2)), center)).^2)); % calculate a eucliean distance for all the 80k vertices at once
                                jammed_vertices_temp = find(distance < sp_kernel);                                                                  % detect the only vertices that are distant within 2mm from the center
                                [a b] = intersect(jammed_vertices_temp, nbr(idx_vertex(iv), :));                                                    % if they are just first order neighbourhoods, then
                                if(~isempty(b))
                                    jammed_vertices_temp(b) = [];                                                                                   % then don't consider as a candidate so exclude them in the cue
                                end
                                jammed_vertices_temp(jammed_vertices_temp == idx_vertex(iv)) = [];                                                  % also exclude the center itself
                                geoDist_temp(jammed_vertices_temp) = min(geoDist.surf(i + 5, [idx_vertex(iv) jammed_vertices_temp], p));
                            end
                            geoDist_temp2{p} = geoDist_temp;
                        end
                        
                        for p = 1 : size(case_num_pat, 1)
                            eval([ 'geoDist_' Modality{m} '_cont_new{' num2str(j) '}.surf(' num2str(i + 5) ', :, ' num2str(p) ') = geoDist_temp2{' num2str(p) '};'  ]);
                        end
                        toc
                        
                    end
                    
                end
                disp(['case : ' case_num_cont{j} ' done']);
                
            end
            
            delete(gcp);
            
        end
                
        for disease_control = 1
            
            OUTPATH          = '/local_raid/seokjun/01_project/04_IntracorticalAnalysis/03_result/97_3T_TLE_HS/';
            fMRI_surfpath    = '/host/weka/export02/data/min/fMRI/DPARSF_processing/DPARSFA_scrub/02_without_mean_sig_regout/bbr_rest/';
                        
            for j = 1 : size(case_num_pat, 1)
                disp(['patient: ', num2str(j)]);
                idx = find(strcmp(case_num_pat, case_num_pat{j}));
                if(exist([ lesion_label_dir '/' Prefix_pat '_' case_num_pat{j} '_label_union_left_rsl.txt' ], 'file'))
                    lesion_label_data_temp = SurfStatReadData( { [ lesion_label_dir '/' Prefix_pat '_' case_num_pat{j} '_label_union_left_rsl.txt' ], ...
                        [ lesion_label_dir '/' Prefix_pat '_' case_num_pat{j} '_label_union_right_rsl.txt' ] } );
                    eval([ 'lesion_label_data(', num2str(j), ', :) = lesion_label_data_temp;']);
                end
            end
            
            parpool(24);
            
            % 1) The original method
            for c = 1 : size(case_num_disease_cont, 1)
                
                for m = 1 : size(Modality, 2)
                    
                    for T1_flair_dti = 1
                        
                        if(strcmp(Modality{m}, 't1') | strcmp(Modality{m}, 'flair') | strcmp(Modality{m}, 'dti'))
                            if(strcmp(Modality{m}, 't1') | strcmp(Modality{m}, 'flair'))
                                eval([ 'geoDist_' Modality{m} '_disease_cont{' num2str(c) '}.surf = zeros(NumIntSurf + NumSubSurf, 81924,  size(case_num_pat, 1));' ]);
                            end
                            if(strcmp(Modality{m}, 'dti'))
                                eval([ 'geoDist_' Modality{m} '_disease_cont{' num2str(c) '}.surf = zeros(NumIntSurf + 2 + NumSubSurf, 81924,  size(case_num_pat, 1));' ]);
                            end
                            
                            surf_ind_cont{c}.surf = cell(NumIntSurf + 2 + NumSubSurf, 1);
                            
                            for i = 1 : NumIntSurf + 2
                                if(i == 1)
                                    basename  = [ Prefix_cont '_' case_num_disease_cont{c} '_gray_surface' ];
                                    temp_str = 'gm';
                                elseif(i == NumIntSurf + 2)
                                    basename  = [ Prefix_cont '_' case_num_disease_cont{c} '_white_surface' ];
                                    temp_str = 'wm';
                                else
                                    basename  = [ Prefix_cont '_' case_num_disease_cont{c} '_intracortical_surface_' num2str(i-1) ];
                                    temp_str = num2str(i-1);
                                end
                                
                                surf_l_name = [ OUTPATH case_num_disease_cont{c} '/surfaces/' basename '_left_'  num2str(NumMesh) postfix_surf '_' Modality{m} '.obj' ];
                                surf_r_name = [ OUTPATH case_num_disease_cont{c} '/surfaces/' basename '_right_' num2str(NumMesh) postfix_surf '_' Modality{m} '.obj' ];
                                surf_l_rsl_name = [ OUTPATH case_num_disease_cont{c} '/surfaces/' basename '_left_'  num2str(NumMesh) postfix_surf '_' Modality{m} '_rsl.obj' ];
                                surf_r_rsl_name = [ OUTPATH case_num_disease_cont{c} '/surfaces/' basename '_right_' num2str(NumMesh) postfix_surf '_' Modality{m} '_rsl.obj' ];
                                xfm_l_name  = [ OUTPATH case_num_disease_cont{c} '/xfm/' Prefix_cont '_' case_num_disease_cont{c} '_left_surfmap.sm' ];
                                xfm_r_name  = [ OUTPATH case_num_disease_cont{c} '/xfm/' Prefix_cont '_' case_num_disease_cont{c} '_right_surfmap.sm' ];
                                
                                if(exist(surf_l_rsl_name, 'file') & exist(surf_r_rsl_name, 'file'))
                                    disp([ basename '_' Modality{m} ' exists.' ]);
                                    surf   = SurfStatReadSurf({ surf_l_rsl_name, surf_r_rsl_name });
                                    surf_ind_cont{c}.surf{i}.tri = surf.tri;
                                    surf_ind_cont{c}.surf{i}.coord = surf.coord;
                                    surf.nbr = nbr;
                                    eval([ 'surf_ind_' Modality{m} '_disease_cont{' num2str(c) '}.surf{' num2str(i) '}    = surf_ind_cont{' num2str(c) '}.surf{' num2str(i) '};' ]);
                                    geoDist_cont_temp = cell(size(case_num_pat, 1), 1);
                                    parfor j = 1 : size(case_num_pat, 1)
                                        disp(['control: ', num2str(c), ', modality: ', num2str(m), ' surface:, ', num2str(i), ', patient: ', num2str(j)]);
                                        geoDist_cont_temp{j} = surfGeoDistHong(surf, lesion_label_data(j, :));
                                    end
                                    geoDist_cont_temp = cell2mat(geoDist_cont_temp);
                                    eval([ 'geoDist_'  Modality{m} '_disease_cont{' num2str(c) '}.surf(' num2str(i) ', :, :) = geoDist_cont_temp''; ' ]);
                                    continue;
                                end
                                
                                if(~exist(surf_l_name, 'file') | ~exist(surf_r_name, 'file'))
                                    continue;
                                end
                                
                                disp([ basename '_' Modality{m} ' is resampled now!' ]);
                                [r, s] = system(['sphere_resample_obj ' surf_l_name ' ' xfm_l_name ' ' surf_l_rsl_name ]);
                                [r, s] = system(['sphere_resample_obj ' surf_r_name ' ' xfm_r_name ' ' surf_r_rsl_name ]);
                                
                                surf   = SurfStatReadSurf({ surf_l_rsl_name, surf_r_rsl_name });
                                surf_ind_cont{c}.surf{i}.tri   = surf.tri;
                                surf_ind_cont{c}.surf{i}.coord = surf.coord;
                                surf.nbr = nbr;
                                eval([ 'surf_ind_' Modality{m} '_disease_cont{' num2str(c) '}.surf{' num2str(i) '}    = surf_ind_cont{' num2str(c) '}.surf{' num2str(i) '};' ]);
                                geoDist_cont_temp = cell(size(case_num_pat, 1), 1);
                                parfor j = 1 : size(case_num_pat, 1)
                                    disp(['control: ', num2str(c), ', modality: ', num2str(m), ' surface:, ', num2str(i), ', patient: ', num2str(j)]);
                                    geoDist_cont_temp{j} = surfGeoDistHong(surf, lesion_label_data(j, :));
                                end
                                geoDist_cont_temp = cell2mat(geoDist_cont_temp);
                                eval([ 'geoDist_'  Modality{m} '_disease_cont{' num2str(c) '}.surf(' num2str(i) ', :, :) = geoDist_cont_temp''; ' ]);
                                
                            end
                            
                            for i = 1 : NumSubSurf
                                basename = [ Prefix_cont '_' case_num_disease_cont{c} '_white_surface_' num2str(i) ];
                                temp_str = num2str(i);
                                
                                surf_l_name = [ OUTPATH case_num_disease_cont{c} '/surfaces/' basename '_left_'  num2str(NumMesh) postfix_surf '_' Modality{m} '.obj' ];
                                surf_r_name = [ OUTPATH case_num_disease_cont{c} '/surfaces/' basename '_right_' num2str(NumMesh) postfix_surf '_' Modality{m} '.obj' ];
                                surf_l_rsl_name = [ OUTPATH case_num_disease_cont{c} '/surfaces/' basename '_left_'  num2str(NumMesh) postfix_surf '_' Modality{m} '_rsl.obj' ];
                                surf_r_rsl_name = [ OUTPATH case_num_disease_cont{c} '/surfaces/' basename '_right_' num2str(NumMesh) postfix_surf '_' Modality{m} '_rsl.obj' ];
                                xfm_l_name  = [ OUTPATH case_num_disease_cont{c} '/xfm/' Prefix_cont '_' case_num_disease_cont{c} '_left_surfmap.sm' ];
                                xfm_r_name  = [ OUTPATH case_num_disease_cont{c} '/xfm/' Prefix_cont '_' case_num_disease_cont{c} '_right_surfmap.sm' ];
                                
                                if(exist(surf_l_rsl_name, 'file') & exist(surf_r_rsl_name, 'file'))
                                    disp([ basename '_' Modality{m} ' exists.' ]);
                                    surf   = SurfStatReadSurf({ surf_l_rsl_name, surf_r_rsl_name });
                                    surf_ind_cont{c}.surf{i+5}.tri = surf.tri;
                                    surf_ind_cont{c}.surf{i+5}.coord = surf.coord;
                                    surf.nbr = nbr;
                                    eval([ 'surf_ind_' Modality{m} '_disease_cont{' num2str(c) '}.surf{' num2str(i+5) '} = surf_ind_cont{' num2str(c) '}.surf{' num2str(i+5) '};' ]);
                                    geoDist_cont_temp = cell(size(case_num_pat, 1), 1);
                                    parfor j = 1 : size(case_num_pat, 1)
                                        disp(['control: ', num2str(c), ', modality: ', num2str(m), ' surface:, ', num2str(i+5), ', patient: ', num2str(j)]);
                                        geoDist_cont_temp{j} = surfGeoDistHong(surf, lesion_label_data(j, :));
                                    end
                                    geoDist_cont_temp = cell2mat(geoDist_cont_temp);
                                    eval([ 'geoDist_'  Modality{m} '_disease_cont{' num2str(c) '}.surf(' num2str(i+5) ', :, :) = geoDist_cont_temp''; ' ]);
                                    continue;
                                end
                                
                                if(~exist(surf_l_name, 'file') | ~exist(surf_r_name, 'file'))
                                    continue;
                                end
                                disp([ basename '_' Modality{m} ' is resampled now!' ]);
                                [r, s] = system(['sphere_resample_obj ' surf_l_name ' ' xfm_l_name ' ' surf_l_rsl_name ]);
                                [r, s] = system(['sphere_resample_obj ' surf_r_name ' ' xfm_r_name ' ' surf_r_rsl_name ]);
                                
                                surf   = SurfStatReadSurf({ surf_l_rsl_name, surf_r_rsl_name });
                                surf_ind_cont{c}.surf{i+5}.tri = surf.tri;
                                surf_ind_cont{c}.surf{i+5}.coord = surf.coord;
                                surf.nbr = nbr;
                                eval([ 'surf_ind_' Modality{m} '_disease_cont{' num2str(c) '}.surf{' num2str(i+5) '} = surf_ind_cont{' num2str(c) '}.surf{' num2str(i+5) '};' ]);
                                geoDist_cont_temp = cell(size(case_num_pat, 1), 1);
                                parfor j = 1 : size(case_num_pat, 1)
                                    disp(['control: ', num2str(c), ', modality: ', num2str(m), ' surface:, ', num2str(i+5), ', patient: ', num2str(j)]);
                                    geoDist_cont_temp{j} = surfGeoDistHong(surf, lesion_label_data(j, :));
                                end
                                geoDist_cont_temp = cell2mat(geoDist_cont_temp);
                                eval([ 'geoDist_'  Modality{m} '_disease_cont{' num2str(c) '}.surf(' num2str(i+5) ', :, :) = geoDist_cont_temp''; ' ]);
                            end
                            
                        end
                    end
                    
                    for fMRI = 1
                        
                        if(strcmp(Modality{m}, 'fMRI'))
                            
                            eval([ 'geoDist_' Modality{m} '_disease_cont{' num2str(c) '}.surf = zeros(size(case_num_pat, 1), 81924);' ]);
                            
                            surf_ind_cont{c}.surf = cell(1, 1);
                            
                            surf_l_name = [ fMRI_surfpath '/' Prefix_cont '_' case_num_disease_cont{c} '/' Prefix_cont '_' case_num_disease_cont{c} '_mid_surface_left_' num2str(NumMesh)  '_rest_bbr.obj' ];
                            surf_r_name = [ fMRI_surfpath '/' Prefix_cont '_' case_num_disease_cont{c} '/' Prefix_cont '_' case_num_disease_cont{c} '_mid_surface_right_' num2str(NumMesh) '_rest_bbr.obj' ];
                            surf_l_rsl_name = [ fMRI_surfpath '/' Prefix_cont '_' case_num_disease_cont{c} '/' Prefix_cont '_' case_num_disease_cont{c} '_mid_surface_left_' num2str(NumMesh)  '_rest_bbr_rsl.obj' ];
                            surf_r_rsl_name = [ fMRI_surfpath '/' Prefix_cont '_' case_num_disease_cont{c} '/' Prefix_cont '_' case_num_disease_cont{c} '_mid_surface_right_' num2str(NumMesh) '_rest_bbr_rsl.obj' ];
                            
                            xfm_l_name  = [ OUTPATH case_num_disease_cont{c} '/xfm/' Prefix_cont '_' case_num_disease_cont{c} '_left_surfmap.sm' ];
                            xfm_r_name  = [ OUTPATH case_num_disease_cont{c} '/xfm/' Prefix_cont '_' case_num_disease_cont{c} '_right_surfmap.sm' ];
                            
                            if(exist(surf_l_rsl_name, 'file') & exist(surf_r_rsl_name, 'file'))
                                disp([ Prefix_cont '_' case_num_disease_cont{c} '_mid_surface_rest_bbr_rsl exists.' ]);
                                surf   = SurfStatReadSurf({ surf_l_rsl_name, surf_r_rsl_name });
                                surf_ind_cont{c}.surf{1}.tri = surf.tri;
                                surf_ind_cont{c}.surf{1}.coord = surf.coord;
                                surf.nbr = nbr;
                                eval([ 'surf_ind_' Modality{m} '_disease_cont{' num2str(c) '}.surf{1} = surf_ind_cont{' num2str(c) '}.surf{1};' ]);
                                geoDist_cont_temp = cell(size(case_num_pat, 1), 1);
                                parfor j = 1 : size(case_num_pat, 1)
                                    disp(['control: ', num2str(c), ', modality: ', num2str(m), ' surface:, mid, patient: ', num2str(j)]);
                                    geoDist_cont_temp{j} = surfGeoDistHong(surf, lesion_label_data(j, :));
                                end
                                geoDist_cont_temp = cell2mat(geoDist_cont_temp);
                                eval([ 'geoDist_'  Modality{m} '_disease_cont{' num2str(c) '}.surf = geoDist_cont_temp; ' ]);
                                continue;
                            end
                            
                            if(~exist(surf_l_name, 'file') | ~exist(surf_r_name, 'file'))
                                continue;
                            end
                            
                            disp([ Prefix_cont '_' case_num_disease_cont{c} '_mid_surface_rest_bbr is resampled now!' ]);
                            [r, s] = system(['/data/noel/noel6/seokjun/02_development/sphere_resample_obj ' surf_l_name ' ' xfm_l_name ' ' surf_l_rsl_name ]);
                            [r, s] = system(['/data/noel/noel6/seokjun/02_development/sphere_resample_obj ' surf_r_name ' ' xfm_r_name ' ' surf_r_rsl_name ]);
%                             [r, s] = system(['sphere_resample_obj ' surf_l_name ' ' xfm_l_name ' ' surf_l_rsl_name ]);
%                             [r, s] = system(['sphere_resample_obj ' surf_r_name ' ' xfm_r_name ' ' surf_r_rsl_name ]);
                            
                            surf   = SurfStatReadSurf({ surf_l_rsl_name, surf_r_rsl_name });
                            surf_ind_cont{c}.surf{1}.tri = surf.tri;
                            surf_ind_cont{c}.surf{1}.coord = surf.coord;
                            surf.nbr = nbr;
                            eval([ 'surf_ind_' Modality{m} '_disease_cont{' num2str(c) '}.surf{1} = surf_ind_cont{' num2str(c) '}.surf{1};' ]);
                            geoDist_cont_temp = cell(size(case_num_pat, 1), 1);
                            parfor j = 1 : size(case_num_pat, 1)
                                disp(['control: ', num2str(c), ', modality: ', num2str(m), ' surface:, mid, patient: ', num2str(j)]);
                                geoDist_cont_temp{j} = surfGeoDistHong(surf, lesion_label_data(j, :));
                            end
                            geoDist_cont_temp = cell2mat(geoDist_cont_temp);
                            eval([ 'geoDist_'  Modality{m} '_disease_cont{' num2str(c) '}.surf = geoDist_cont_temp; ' ]);
                            
                        end
                        
                    end
                    
                end
                
                disp(['case : ' case_num_disease_cont{c} ' done']);
                
            end
            
            delete(gcp);
            
            % 2) Modified geodesic distance
            % Problem: An originally calculated geodesic distance on subcortical white matter surface has an intrinsic limitation
            % at the narrow white matter region where two surfaces are jammed together very closely.
            % It's because although the original geodesic distance is gradually growing along the surface as far from the lesion,
            % sampling the value is being done still in same area.
            % Solution: 1) We apply a small sphere kernel at every vertex to detect such jammed surface areas.
            %              If more than two vertices are included in the sphere, we assume the surfaces at that location are likely jammed each other.
            %           2) We assign minimal geodesic distance to all detected vertices.
            sp_kernel = 2;
            surf = surfGetNeighborsHong(surf);
            nbr = surf.nbr;
            load('individual_surface_multimodal_disease_cont.mat');
            
            parpool(24);
            
            geoDist_t1_disease_cont_new = geoDist_t1_disease_cont;
            geoDist_flair_disease_cont_new = geoDist_flair_disease_cont;
            geoDist_dti_disease_cont_new = geoDist_dti_disease_cont;
            geoDist_fMRI_disease_cont_new = geoDist_fMRI_disease_cont;
            for j = 1 : size(case_num_disease_cont, 1)
                
                for m = 1 : size(Modality, 2)
                    
                    eval([ 'geoDist  = geoDist_' Modality{m} '_disease_cont{' num2str(j) '};'  ]);
                    eval([ 'surf_ind = surf_ind_' Modality{m} '_disease_cont{' num2str(j) '};'  ]);
                    
                    if(strcmp(Modality{m}, 't1') | strcmp(Modality{m}, 'flair'))
                        NumSubSurf_temp = 3;
                    elseif(strcmp(Modality{m}, 'dti'))
                        NumSubSurf_temp = 5;
                    elseif(strcmp(Modality{m}, 'fMRI'))
                        NumSubSurf_temp = 0;
                    end
                    
                    for i = 1 : NumSubSurf_temp
                        
                        tic
                        geoDist_temp2 = cell(size(case_num_pat, 1), 1);
                        parfor p = 1 : size(case_num_pat, 1)
                            disp(['cont : ', case_num_disease_cont{j}, ', mod : ', Modality{m}, ', surf : ', num2str(i+5) ', pat : ', case_num_pat{p}]);
                            idx_vertex = find(geoDist.surf(i + 5, :, p) < 25 & geoDist.surf(i + 5, :, p) > 0);                                      % limited the circal range where we detect a jammed surface area within 50mm
                            jammed_vertices = [];
                            geoDist_temp = geoDist.surf(i + 5, :, p);
                            for iv = 1 : size(idx_vertex, 2)
                                center = surf_ind.surf{i + 5}.coord(:, idx_vertex(iv));                                                             % for those vertices
                                distance = sqrt(sum((surf_ind.surf{i + 5}.coord - kron(ones(1, size(surf_ind.surf{i + 5}.coord, 2)), center)).^2)); % calculate a eucliean distance for all the 80k vertices at once
                                jammed_vertices_temp = find(distance < sp_kernel);                                                                  % detect the only vertices that are distant within 2mm from the center
                                [a b] = intersect(jammed_vertices_temp, nbr(idx_vertex(iv), :));                                                    % if they are just first order neighbourhoods, then
                                if(~isempty(b))
                                    jammed_vertices_temp(b) = [];                                                                                   % then don't consider as a candidate so exclude them in the cue
                                end
                                jammed_vertices_temp(jammed_vertices_temp == idx_vertex(iv)) = [];                                                  % also exclude the center itself
                                geoDist_temp(jammed_vertices_temp) = min(geoDist.surf(i + 5, [idx_vertex(iv) jammed_vertices_temp], p));
                            end
                            geoDist_temp2{p} = geoDist_temp;
                        end
                        
                        for p = 1 : size(case_num_pat, 1)
                            eval([ 'geoDist_' Modality{m} '_disease_cont_new{' num2str(j) '}.surf(' num2str(i + 5) ', :, ' num2str(p) ') = geoDist_temp2{' num2str(p) '};'  ]);
                        end
                        toc
                        
                    end
                    
                end
                disp(['case : ' case_num_disease_cont{j} ' done']);
                
            end
            
            delete(gcp);
            
        end
        
    end
    
end
