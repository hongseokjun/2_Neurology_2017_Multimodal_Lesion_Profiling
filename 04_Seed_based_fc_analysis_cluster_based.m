clear all
close all

%% path configuration
addpath('/data/noel/noel2/local/matlab/surfstat_chicago');
addpath(genpath('/data/noel/noel6/seokjun/00_local_raid/03_downloads/tSNE_matlab/'));

%% configure analysis parameters
meansig_cov = 'wbsignocov'; % wbsigcov or wbsignocov
smooth_method = 'minc';   % minc or surfstat
printfigs = 0;            % 1 (print) or 0 (no print)
load_mat = 1;             % 1 (load pre-saved files) or 0 (newly read all files)
save_file = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

%% Analyses configuration is being setting-up ...
%% 1) setup directories
for setup_directories = 1
    if(strcmp(meansig_cov, 'wbsigcov'))
        WDIR = '/host/gypsy/local_raid/seokjun/01_project/05_NoelRest/03_Result/01_DPARSF_processing/DPARSFA_scrub/01_with_mean_sig_regout/';
    else
        WDIR = '/host/gypsy/local_raid/seokjun/01_project/05_NoelRest/03_Result/01_DPARSF_processing/DPARSFA_scrub/02_without_mean_sig_regout/';
    end
    SDIR = '/local_raid/seokjun/01_project/05_NoelRest/01_Analysis/';
    RDIR = '/host/gypsy/local_raid/seokjun/01_project/05_NoelRest/03_Result/03_functional_network_topology/';
    CDIR = '/local_raid/seokjun/01_project/05_NoelRest/00_CaseList/';
    LDIR = '/data/noel/noel6/CTFCD-1.2.0_64/Lesion/lesion_surf/';
    PDIR = '/local_raid/seokjun/01_project/05_NoelRest/01_Analysis/parcellation/';
    
    addpath(SDIR);
end

%% 2) read surfaces
for readsurf = 1
    
    surf  = SurfStatAvSurf({[ SDIR '/' 'surf_reg_model_left.obj' ],[ SDIR '/' 'surf_reg_model_right.obj' ]});
    mask  = SurfStatMaskCut(surf);
    
    if(printfigs)
        f     = figure; BoSurfStatViewData(mask,surf,'');
    end
    
    isurf = SurfStatInflate(surf,0.3);
    Curv  = SurfStatReadData({[ SDIR '/' 'surf_reg_model_left.txt' ] ,[ SDIR '/' 'surf_reg_model_left.txt' ]});
    
    if(printfigs)
        f     = figure; BoSurfStatViewData(sign(Curv),isurf,'');
        colormap([0.6 0.6 0.6;0.8 0.8 0.8]);
    end
    
    rois_L = SurfStatReadData1([ PDIR '/surf_reg_model_left_parcellation_kmeans.txt' ]);
    
    %% Left: 0-583 | Right: 584-1167
    rois = [ rois_L, rois_L+584];
    rois(rois==584) = 0;
    rois(rois>584) = rois(rois>584) - 1;    
    
    if(printfigs)
        f     = figure; BoSurfStatViewData(rois,isurf,'');
        colormap([jet;jet])
    end
    surf_lowres = SurfStatAvSurf({ 'surf_reg_model_left_20480.obj', 'surf_reg_model_right_20480.obj' });
end

%% 3) setup colormap
for colormaps = 1
    
      blackblue = [ zeros(1,3)*0.8; 
                    zeros(127,1)   (0:126)'/127   ones(127,1)];
      blackblue = flipud(blackblue);

      blue      = [ ones(1,3)*0.8; 
                    zeros(127,1)   (0:126)'/127   ones(127,1)];
      blue      = flipud(blue);

      red       = [  ones(1,3)*0.8; ...
                     ones(64,1) linspace(0,253,64)'/254 zeros(64,1);...
                     ones(64,1) ones(64,1) linspace(0,253,64)'/254];         
                   
      mc2       = [ makeColorMap([1 1 0],[1 0 0],[0 0 0],512/4); ...
                    makeColorMap([0 0 0],[0.7 0.7 0.7],[0 0 0],812/4); ...
                    makeColorMap([0 0 0],[0 0 1],[0 1 1],512/4)];
      mc2       = flipud(mc2); 
      
end

%% 4) read csv
for readcsv = 1
    
    fid   = fopen([ CDIR 'data_matlab_func_new.csv' ]) ;
    C     = textscan(fid,'%s%f%s%s%s%n','delimiter',',','CollectOutput',1);    
    CODE  = C{1}(:, 1); 
    AGE   = C{2}(:, 1);
    SEX   = C{3}(:, 1);
    GROUP = C{3}(:, 2); 
    Histo = C{3}(strcmp(GROUP, 'FCD'), 3);    

    PREFIX                      = GROUP; 
    PREFIX(strcmp(GROUP,'NC'))  = {'TLE'}; 
    PREFIX(strcmp(GROUP,'FCD')) = {'mcd'}; 
    FOLDER                      = GROUP; 
    FOLDER(strcmp(GROUP,'NC'))  = {WDIR}; 
    FOLDER(strcmp(GROUP,'FCD')) = {WDIR}; 
    
    % Age and sex test
    % 1-1) FCD vs. NC: AGE
    [h,p,ci,stats] = ttest2(AGE(strcmp(GROUP, 'FCD')), AGE(strcmp(GROUP, 'NC')), 0.05, 'both')
    [ mean(AGE(strcmp(GROUP, 'FCD'))) std(AGE(strcmp(GROUP, 'FCD'))) ;
      mean(AGE(strcmp(GROUP, 'NC')))  std(AGE(strcmp(GROUP, 'NC'))) ]
  
    % 1-2) FCD vs. NC: SEX
    [table,chi2,p,labels] = crosstab([sum(strcmp(SEX, 'f')&strcmp(GROUP, 'FCD')); sum(strcmp(SEX, 'm')&strcmp(GROUP, 'FCD'))], ...
                                     [sum(strcmp(SEX, 'f')&strcmp(GROUP, 'NC')); sum(strcmp(SEX, 'm')&strcmp(GROUP, 'NC'))])
          
    % 2-1) FCD subtupe vs. NC: AGE
    [h,p,ci,stats] = ttest2(AGE(strcmp(Histo, 'IIA')), AGE(strcmp(GROUP, 'NC')), 0.05, 'both')
    [ mean(AGE(strcmp(Histo, 'IIA'))) std(AGE(strcmp(Histo, 'IIA'))) ;
      mean(AGE(strcmp(GROUP, 'NC')))  std(AGE(strcmp(GROUP, 'NC'))) ]
  
    [h,p,ci,stats] = ttest2(AGE(strcmp(Histo, 'IIB')), AGE(strcmp(GROUP, 'NC')), 0.05, 'both')
    [ mean(AGE(strcmp(Histo, 'IIB'))) std(AGE(strcmp(Histo, 'IIB'))) ;
      mean(AGE(strcmp(GROUP, 'NC')))  std(AGE(strcmp(GROUP, 'NC'))) ]  
  
    % 2-2) FCD subtupe vs. NC: SEX
    [table,chi2,p,labels] = crosstab([sum(strcmp(SEX, 'f')&strcmp(C{3}(:, 3), 'IIA')); sum(strcmp(SEX, 'm')&strcmp(C{3}(:, 3), 'IIA'))], ...
                                     [sum(strcmp(SEX, 'f')&strcmp(GROUP, 'NC')); sum(strcmp(SEX, 'm')&strcmp(GROUP, 'NC'))])   
    [table,chi2,p,labels] = crosstab([sum(strcmp(SEX, 'f')&strcmp(C{3}(:, 3), 'IIB')); sum(strcmp(SEX, 'm')&strcmp(C{3}(:, 3), 'IIB'))], ...
                                     [sum(strcmp(SEX, 'f')&strcmp(GROUP, 'NC')); sum(strcmp(SEX, 'm')&strcmp(GROUP, 'NC'))])
                                 
    num_of_cont    = sum(strcmp(GROUP, 'NC'))
    num_of_FCD     = sum(strcmp(GROUP, 'FCD'))
    num_of_FCDIIA  = sum(strcmp(Histo, 'IIA'))
    num_of_FCDIIB  = sum(strcmp(Histo, 'IIB'))
    
end

%% 5) downsample surface and time series
data_temp = zeros(1, 81924);
downsample_num = 2562;
data_temp([ 1:downsample_num 40963:40962+downsample_num ]) = 1;
figure; SurfStatView1(data_temp, surf); axis off; cameramenu; material dull;
view(-90, 0); light('position', [-1.0872e+03 -256.5311 183.7751]);

ts_set = zeros(145, 81924, num_of_cont+num_of_FCD);
for i = 1 : length(GROUP)
    ts_name   = [FOLDER{i} '/bbr_rest/' PREFIX{i} '_' CODE{i} '/surf_ts/' PREFIX{i} '_' CODE{i} '_4D_sm_5_rsl.mgh'];
    ts_set(:, :, i) = SurfStatReadData_mgh1(ts_name);
end

order = 1;
edg=SurfStatEdg(surf);
nbr_L = cell(downsample_num, 1);
parfor i = 1 : downsample_num
    i
    nbr_L(i) = { find_vertex_multineighbors_by_order(surf, i, order, edg)' };
end

nbr_R = cell(downsample_num, 1);
parfor i = 1 : downsample_num
    i
    nbr_R(i) = { find_vertex_multineighbors_by_order(surf, i+size(surf.coord, 2)/2, order, edg)' };
end

nbr = [ nbr_L; nbr_R ];

ts_set_resampled = zeros(size(ts_set, 1), downsample_num*2, length(GROUP));
for i = 1 : length(GROUP)
    
    i
   for j = 1 : downsample_num*2
        ts_set_resampled(:, j, i) = mean(ts_set(:, nbr{j}, i), 2);
    end

end

surf_temp = SurfStatReadSurf({[ '/local_raid/seokjun/01_project/05_NoelRest/99_temp/ellipsoid_' num2str((downsample_num-2)*2) '.obj' ], ...
                              [ '/local_raid/seokjun/01_project/05_NoelRest/99_temp/ellipsoid_' num2str((downsample_num-2)*2) '.obj' ] });
surf_lowres.tri = surf_temp.tri;
surf_lowres.coord = surf.coord(:, [ 1:downsample_num 40963:40963+downsample_num-1]);
mask_lowres = SurfStatMaskCut(surf_lowres);

%% 6) compute a fc map of controls
conn_control_fc_matrices = cell(num_of_FCD, 1);
for i = 1 : num_of_FCD

    %%  0) lesion label setting
    disp(['case ' CODE{i} ' start!']);
    if(exist([ LDIR '/' PREFIX{i} '_' CODE{i} '_label_union_left_rsl.txt' ], 'file'))
        lesion_label_data = SurfStatReadData( { ...
            [ LDIR '/' PREFIX{i} '_' CODE{i} '_label_union_left_rsl.txt'], ...
            [ LDIR '/' PREFIX{i} '_' CODE{i} '_label_union_right_rsl.txt' ] } );
    end
    lesion_label_data = lesion_label_data([1:downsample_num 40963:40963+downsample_num-1]);
    lesion_label_data_set(i, :) = lesion_label_data;        
    lesion_curr = find(lesion_label_data_set(i, :)~=0);
    
    %%  1) At each control, arcoss every patient, within same anatomical area to the lesion, compute vertex-wise functional connectivity to the rest of brain -> generate a fc matrix at each control
    conn_lesional_fc_matrices = zeros(sum(lesion_label_data_set(i, :)~=0), downsample_num*2, num_of_cont);
    for j = 1 : num_of_cont
        
        conn_lesional_fc_matrices(:, :, j) = transpose(corr(ts_set_resampled(:, :, num_of_FCD+j), ts_set_resampled(:, lesion_curr, num_of_FCD+j)));
        
    end
    
    %%  2) Average these matrices across controls -> generate one matrix for average pattern of normal fc regard to the lesion
    average_conn_lesional_fc_matrix = mean(conn_lesional_fc_matrices, 3);
    average_conn_lesional_fc_matrix(isnan(average_conn_lesional_fc_matrix)) = 0;
    
    %%  3) Visually check the average matrix to already detect grouped patterns
    
    %%  4) Based on a group-level connectivity matrix, cluster the lesional vertices into  distinctive patterns using k-mean clustering
    
    %%  5) Decide the optimized number of clusters and subdivide the lesional vertices
    %      and recalculate seed-based connectivity based on subdivided clusters
    if(sum(lesion_label_data_set(i, :)~=0)>10)    
        
        eva = evalclusters(average_conn_lesional_fc_matrix,'kmeans','Silhouette','KList',[1:10]);
        clusters = eva.OptimalY;
        
        conn_lesional_fc_matrices_temp = zeros(num_of_cont, size(surf_lowres.coord, 2), eva.OptimalK);
        subdivied_lesion = zeros(1, size(surf_lowres.coord, 2));
        for j = 1 : eva.OptimalK
            
            subdivied_lesion(lesion_curr(clusters==j)) = j;
            clus_id = j;
            
            for k = 1 : num_of_cont
                
                conn_lesional_fc_matrices_temp(k, :, j) = transpose(corr(ts_set_resampled(:, :, num_of_FCD+k), mean(ts_set_resampled(:, subdivied_lesion==clus_id, num_of_FCD+k), 2)));
                
            end
            
        end
        
        conn_control_fc_matrices(i) = { conn_lesional_fc_matrices_temp };
        
    else
        
        conn_lesional_fc_matrices_temp = zeros(num_of_cont, size(surf_lowres.coord, 2), 1);
        for k = 1 : num_of_cont
                
                conn_lesional_fc_matrices_temp(k, :, 1) = transpose(corr(ts_set_resampled(:, :, num_of_FCD+k), mean(ts_set_resampled(:, lesion_label_data_set(i, :)~=0, num_of_FCD+k), 2)));
                
        end
        
        conn_control_fc_matrices(i) = { conn_lesional_fc_matrices };
        
    end
    
end