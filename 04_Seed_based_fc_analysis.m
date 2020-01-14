clear all
close all
addpath('/data/noel/noel2/local/matlab/surfstat_chicago');

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

%% 5) read network data
for read_network_data = 1
    
    % load network data
    for loaddata = 1
        
        ts_set = zeros(145, 81924, num_of_cont+num_of_FCD);
        for i = 1 : length(GROUP)
            ts_name   = [FOLDER{i} '/bbr_rest/' PREFIX{i} '_' CODE{i} '/surf_ts/' PREFIX{i} '_' CODE{i} '_4D_sm_5_rsl.mgh'];
            ts_set(:, :, i) = SurfStatReadData_mgh1(ts_name);
        end
        
        %% 1) change the data precision format to single from double (to spare the memory)
        %% 2) downsample it (to 2562 or 10242)
        downsample_vert = 2562;
        surf_temp = SurfStatReadSurf({[ '/local_raid/seokjun/01_project/05_NoelRest/99_temp/ellipsoid_' num2str((downsample_vert-2)*2) '.obj' ], [ '/local_raid/seokjun/01_project/05_NoelRest/99_temp/ellipsoid_' num2str((downsample_vert-2)*2) '.obj' ] });
        surf_lowres.tri = surf_temp.tri;
        surf_lowres.coord = surf.coord(:, [ 1:downsample_vert 40963:40963+downsample_vert-1]);
        mask_lowres = SurfStatMaskCut(surf_lowres);
        
        ts_set = single(ts_set);
        ts_set = ts_set(:, [ 1:downsample_vert 40963:40963+downsample_vert-1], :);
                    
        for method_1 = 1
            
            lesion_conn_set                 = cell(num_of_FCD, 1);
            cont_conn_set                   = cell(num_of_FCD, 1);
            
            z_lesion_finalset               = zeros(num_of_FCD, 1);
            internal_z_lesion_finalset      = zeros(num_of_FCD, 1);
            average_z_lesion_conn           = zeros(num_of_FCD, downsample_vert*2);
            internal_z_lesion_conn          = zeros(num_of_FCD, downsample_vert*2);
            average_mean_cont_conn          = zeros(num_of_FCD, downsample_vert*2);
            r_fc_cont_conn_finalset         = zeros(num_of_FCD, num_of_cont);
            
            z_flesion_finalset               = zeros(num_of_FCD, 1);
            internal_z_flesion_finalset      = zeros(num_of_FCD, 1);
            average_z_flesion_conn           = zeros(num_of_FCD, downsample_vert*2);
            internal_z_flesion_conn          = zeros(num_of_FCD, downsample_vert*2);
            average_mean_fcont_conn          = zeros(num_of_FCD, downsample_vert*2);
            r_fc_fcont_conn_finalset          = zeros(num_of_FCD, num_of_cont);
            
            for i = 1:num_of_FCD
                
                disp(['case ' CODE{i} ' start!']);
                if(exist([ LDIR '/' PREFIX{i} '_' CODE{i} '_label_union_left_rsl.txt' ], 'file'))
                    lesion_label_data = SurfStatReadData( { ...
                        [ LDIR '/' PREFIX{i} '_' CODE{i} '_label_union_left_rsl.txt'], ...
                        [ LDIR '/' PREFIX{i} '_' CODE{i} '_label_union_right_rsl.txt' ] } );
                end
                lesion_label_data = lesion_label_data([1:downsample_vert 40963:40963+downsample_vert-1]);
                
                lesion_conn = corr(ts_set(:, lesion_label_data~=0, i), ts_set(:, :, i));
                flesion_conn = 0.5*log((1+lesion_conn)./(1-lesion_conn));
                cont_conn     = single(zeros(sum(lesion_label_data~=0), size(ts_set, 2), num_of_cont));
                fcont_conn    = single(zeros(sum(lesion_label_data~=0), size(ts_set, 2), num_of_cont));
                for j = 1 : num_of_cont
                    fprintf([ num2str(j) ' ']);
                    cont_conn(:, :, j)  = corr(ts_set(:, lesion_label_data~=0, num_of_FCD+j), ts_set(:, :, num_of_FCD+j));
                    fcont_conn(:, :, j) = 0.5*log((1+cont_conn(:, :, j))./(1-cont_conn(:, :, j)));
                end
                fprintf('\n');
                lesion_conn_set{i} = lesion_conn;
                cont_conn_set{i}   = cont_conn;
                
                % r
                mean_cont_conn = mean(cont_conn, 3);
                std_cont_conn  = std(cont_conn, 0, 3);
                z_lesion_conn = (lesion_conn - mean_cont_conn)./std_cont_conn;
                average_z_lesion_conn(i, :) = mean(z_lesion_conn, 1);
                average_z_lesion_conn(i, isnan(average_z_lesion_conn(i, :))) = 0;
                average_z_lesion_conn(i, isinf(average_z_lesion_conn(i, :))) = 0;
                internal_z_lesion_conn(i, :) = (average_z_lesion_conn(i, :) - mean(average_z_lesion_conn(i, :), 2))/std(average_z_lesion_conn(i, :), 0, 2);
                average_mean_cont_conn(i, :) = mean(mean_cont_conn, 1); average_mean_cont_conn(i, isnan(average_mean_cont_conn(i, :))) = 0;
                
                % fishers' z-score
                mean_fcont_conn = mean(fcont_conn, 3);
                std_fcont_conn  = std(fcont_conn, 0, 3);
                z_flesion_conn = (flesion_conn - mean_fcont_conn)./std_fcont_conn;
                average_z_flesion_conn(i, :) = mean(z_flesion_conn, 1);
                average_z_flesion_conn(i, isnan(average_z_flesion_conn(i, :))) = 0;
                average_z_flesion_conn(i, isinf(average_z_flesion_conn(i, :))) = 0;
                internal_z_flesion_conn(i, :) = (average_z_flesion_conn(i, :) - mean(average_z_flesion_conn(i, :), 2))/std(average_z_flesion_conn(i, :), 0, 2);
                average_mean_fcont_conn(i, :) = mean(mean_fcont_conn, 1); average_mean_fcont_conn(i, isnan(average_mean_fcont_conn(i, :))) = 0;
            end
            
            count = 1;
            r_thres_set = 0.20:0.05:0.45;
            
            FCDTypeIIB_z_lesion_finalset = [];
            FCDTypeIIA_z_lesion_finalset = [];
            FCDTypeIIB_z_flesion_finalset = [];
            FCDTypeIIA_z_flesion_finalset = [];
            for r_thres = r_thres_set
                count
                for i = 1:num_of_FCD
                    
                    ROI_cont = average_mean_cont_conn(i, :) >= r_thres;
                    if(sum(ROI_cont) == 0)
                        ROI_cont = average_mean_cont_conn(i, :) >= (max(average_mean_cont_conn(i, :))-0.05);
                    end
                    z_lesion_finalset(i)          = mean(average_z_lesion_conn(i, ROI_cont==1));
                    internal_z_lesion_finalset(i) = mean(internal_z_lesion_conn(i, ROI_cont==1));
                    for j = 1:num_of_cont
                        r_fc_cont_conn_finalset(i, j) = mean(mean(cont_conn_set{i}(:, ROI_cont, j), 1), 2);
                    end
                    
                    ROI_cont = average_mean_fcont_conn(i, :) >= r_thres;
                    z_flesion_finalset(i)          = mean(average_z_flesion_conn(i, ROI_cont==1));
                    internal_z_flesion_finalset(i) = mean(internal_z_flesion_conn(i, ROI_cont==1));
                end
                
                abs_z_lesion_finalset = abs(z_lesion_finalset);
                TYPE='IIB'; FCDTypeIIB_z_lesion_finalset(count, :) = [ mean(abs_z_lesion_finalset(strcmp(Histo, TYPE))) std(abs_z_lesion_finalset(strcmp(Histo, TYPE))) ];
                TYPE='IIA'; FCDTypeIIA_z_lesion_finalset(count, :) = [ mean(abs_z_lesion_finalset(strcmp(Histo, TYPE))) std(abs_z_lesion_finalset(strcmp(Histo, TYPE))) ];
                
                abs_z_flesion_finalset = abs(z_flesion_finalset);
                TYPE='IIB'; FCDTypeIIB_z_flesion_finalset(count, :) = [ mean(abs_z_flesion_finalset(strcmp(Histo, TYPE))) std(abs_z_flesion_finalset(strcmp(Histo, TYPE))) ];
                TYPE='IIA'; FCDTypeIIA_z_flesion_finalset(count, :) = [ mean(abs_z_flesion_finalset(strcmp(Histo, TYPE))) std(abs_z_flesion_finalset(strcmp(Histo, TYPE))) ];
                count = count + 1;
            end
            
            figure; errorbar(((1:length(r_thres_set))+0.2)', flip(FCDTypeIIB_z_lesion_finalset(:, 1), 1), flip(FCDTypeIIB_z_lesion_finalset(:, 2), 1), 'k');
            hold on; errorbar(((1:length(r_thres_set))-0.2)', flip(FCDTypeIIA_z_lesion_finalset(:, 1), 1), flip(FCDTypeIIA_z_lesion_finalset(:, 2), 1), 'r');
            xlabel = num2cell([flip(r_thres_set, 2)]);
            xlabel = [ cell(size(1,1)) xlabel cell(size(1,1)) ];
            set(gca, 'XTickLabel', xlabel); ylim([-0.5 1.5]);
            
            f1=figure; SurfStatView(mean(z_lesion_conn, 1), surf_lowres); SurfStatColLim([-2 2]);
            f2=figure; SurfStatView(mean(lesion_conn, 1), surf_lowres); SurfStatColLim([-1 1]);
            f3=figure; SurfStatView(average_r_conn, surf_lowres); SurfStatColLim([-1 1]);
            f3=figure; SurfStatView(ROI_cont*2, surf_lowres); SurfStatColLim([-1 1]);
            
            f1=figure; SurfStatView(mean(z_fishers_lesion_conn, 1), surf_lowres); SurfStatColLim([-2 2]);
            f2=figure; SurfStatView(mean(fishers_trans_lesion_conn, 1), surf_lowres); SurfStatColLim([-1 1]);
            f3=figure; SurfStatView(mean(fishers_mean_conn, 1), surf_lowres); SurfStatColLim([-1 1]);
            
        end
        for method_2 = 1
            
            %% data read and compute the connectivity
            lesion_conn_set                 = cell(num_of_FCD, 1);
            cont_conn_set                   = cell(num_of_FCD, 1);                        
            
            average_z_lesion_conn           = zeros(num_of_FCD, downsample_vert*2);
            internal_z_lesion_conn          = zeros(num_of_FCD, downsample_vert*2);
            average_mean_cont_conn          = zeros(num_of_FCD, downsample_vert*2);            
            
            average_z_flesion_conn          = zeros(num_of_FCD, downsample_vert*2);
            internal_z_flesion_conn         = zeros(num_of_FCD, downsample_vert*2);
            average_mean_fcont_conn         = zeros(num_of_FCD, downsample_vert*2);
            
            z_cont_conn_set                 = zeros(num_of_cont, downsample_vert*2, num_of_FCD);
            z_fcont_conn_set                = zeros(num_of_cont, downsample_vert*2, num_of_FCD);
            lesion_label_data_set           = zeros(num_of_FCD, downsample_vert*2);
            
            for i = 1 : num_of_FCD
                
                disp(['case ' CODE{i} ' start!']);
                if(exist([ LDIR '/' PREFIX{i} '_' CODE{i} '_label_union_left_rsl.txt' ], 'file'))
                    lesion_label_data = SurfStatReadData( { ...
                        [ LDIR '/' PREFIX{i} '_' CODE{i} '_label_union_left_rsl.txt'], ...
                        [ LDIR '/' PREFIX{i} '_' CODE{i} '_label_union_right_rsl.txt' ] } );
                end
                lesion_label_data = lesion_label_data([1:downsample_vert 40963:40963+downsample_vert-1]);
                lesion_label_data_set(i, :) = lesion_label_data;
                
                %% 1) average time series within the lesion at this given patient
                %% 2) compute the connectivity of the lesion with the rest of the brain
                %% 3) transform the connectivity to Fishers' z-score                
                lesion_ts = mean(ts_set(:, lesion_label_data ~= 0, i), 2);
                lesion_conn = corr(lesion_ts, ts_set(:, :, i));
                flesion_conn = 0.5*(log((1+lesion_conn)./(1-lesion_conn)));
                
                %% 1) average time series within the lesional area at each control
                %% 2) compute the connectivity of the area corresponding to the lesion with the rest of the brain
                %% 3) transform the connectivity to Fishers' z-score
                cont_conn = single(zeros(num_of_cont, length(surf_lowres.coord)));
                fcont_conn = single(zeros(num_of_cont, length(surf_lowres.coord)));
                for j = 1 : num_of_cont
                    fprintf([ num2str(j) ' ']);
                    cont_lesion_ts = mean(ts_set(:, lesion_label_data ~= 0, num_of_FCD+j), 2);
                    cont_conn(j, :) = corr(cont_lesion_ts, ts_set(:, :, num_of_FCD+j));
                    fcont_conn(j, :) = 0.5*(log((1+cont_conn(j, :))./(1-cont_conn(j, :))));                    
                end
                fprintf('\n');
                
                %% save these connecitivity measures
                lesion_conn_set{i} = lesion_conn;
                cont_conn_set{i}   = cont_conn;
  
                %% 1) compute the mean and sd of the connectivity matrices across all controls
                %% 2) compute the z-score of the the connectivity matrices in this patient
                %% 3) 'average_z_lesion_conn' saves the 'z_lesion_conn' across iterations
                %%    'average_mean_cont_conn' saves the 'mean_cont_conn' across iterations
                %% 4) compute internal z-score within this case
                %% ** we do above procedures for both r and Fishers' score
                
                %% r
                mean_cont_conn = mean(cont_conn, 1);
                std_cont_conn  = std(cont_conn, 0, 1);
                z_cont_conn_set(:, :, i) = (cont_conn - repmat(mean_cont_conn, num_of_cont, 1))./repmat(std_cont_conn, num_of_cont, 1);
                z_lesion_conn = (lesion_conn - mean_cont_conn)./std_cont_conn;
                average_z_lesion_conn(i, :) = mean(z_lesion_conn, 1);
                average_z_lesion_conn(i, isnan(average_z_lesion_conn(i, :))) = 0;
                average_z_lesion_conn(i, isinf(average_z_lesion_conn(i, :))) = 0;                
                average_mean_cont_conn(i, :) = mean(mean_cont_conn, 1); average_mean_cont_conn(i, isnan(average_mean_cont_conn(i, :))) = 0;
                internal_z_lesion_conn(i, :) = (average_z_lesion_conn(i, :) - mean(average_z_lesion_conn(i, :), 2))/std(average_z_lesion_conn(i, :), 0, 2);
                
                %% fishers' score
                mean_fcont_conn = mean(fcont_conn, 1);
                std_fcont_conn  = std(fcont_conn, 0, 1);
                z_fcont_conn_set(:, :, i) = (fcont_conn - repmat(mean_fcont_conn, num_of_cont, 1))./repmat(std_fcont_conn, num_of_cont, 1);
                z_flesion_conn = (flesion_conn - mean_fcont_conn)./std_fcont_conn;
                average_z_flesion_conn(i, :) = mean(z_flesion_conn, 1);
                average_z_flesion_conn(i, isnan(average_z_flesion_conn(i, :))) = 0;
                average_z_flesion_conn(i, isinf(average_z_flesion_conn(i, :))) = 0;
                internal_z_flesion_conn(i, :) = (average_z_flesion_conn(i, :) - mean(average_z_flesion_conn(i, :), 2))/std(average_z_flesion_conn(i, :), 0, 2);
                average_mean_fcont_conn(i, :) = mean(mean_fcont_conn, 1); average_mean_fcont_conn(i, isnan(average_mean_fcont_conn(i, :))) = 0;
                
            end
            
            %% statiatical test and visualize the results
            for statistical_test = 1
                
                count = 1;
                r_thres_set = 0.20:0.01:0.70;
                
                r_fc_cont_conn_finalset       = zeros(num_of_FCD, num_of_cont, length(r_thres_set));
                z_lesion_finalset             = zeros(num_of_FCD, length(r_thres_set));
                z_cont_finalset             = zeros(num_of_FCD, num_of_cont, length(r_thres_set));
                internal_z_lesion_finalset    = zeros(num_of_FCD, length(r_thres_set));
                z_flesion_finalset            = zeros(num_of_FCD, length(r_thres_set));
                internal_z_flesion_finalset   = zeros(num_of_FCD, length(r_thres_set));
                z_fcont_finalset             = zeros(num_of_FCD, num_of_cont, length(r_thres_set));
                
                FCDTypeIIB_z_lesion_finalset  = [];
                FCDTypeIIA_z_lesion_finalset  = [];
                FCDTypeIIB_z_flesion_finalset = [];
                FCDTypeIIA_z_flesion_finalset = [];
                
                for r_thres = r_thres_set
                    count                  
                    
                    for i = 1:num_of_FCD
                        
                        ROI_cont = average_mean_cont_conn(i, :) >= r_thres;
                        if(sum(ROI_cont) == 0)
                            disp('hey!');
                            ROI_cont = average_mean_cont_conn(i, :) >= (max(average_mean_cont_conn(i, :))-0.05);
                        end
                        ROI_cont = ROI_cont.*(lesion_label_data_set(i, :)==0);
                        
                        z_lesion_finalset(i, count)          = mean(average_z_lesion_conn(i, ROI_cont==1));
                        internal_z_lesion_finalset(i, count) = mean(internal_z_lesion_conn(i, ROI_cont==1));
                        for j = 1:num_of_cont
                            r_fc_cont_conn_finalset(i, j, count) = mean(mean(cont_conn_set{i}(j, ROI_cont==1), 1), 2);
                            z_cont_finalset(i, j, count) = mean(z_cont_conn_set(j, ROI_cont==1, i), 2);
                        end
                        
                        ROI_cont = average_mean_fcont_conn(i, :) >= r_thres;
                        ROI_cont = ROI_cont.*(lesion_label_data_set(i, :)==0);
                        z_flesion_finalset(i, count)          = mean(average_z_flesion_conn(i, ROI_cont==1));
                        internal_z_flesion_finalset(i, count) = mean(internal_z_flesion_conn(i, ROI_cont==1));
                        for j = 1:num_of_cont                            
                            z_fcont_finalset(i, j, count) = mean(z_fcont_conn_set(j, ROI_cont==1, i), 2);
                        end                        
                        
                    end
                    
                    abs_z_lesion_finalset = abs(z_lesion_finalset(:, count));
                    TYPE='IIB'; FCDTypeIIB_z_lesion_finalset(count, :) = [ mean(abs_z_lesion_finalset(strcmp(Histo, TYPE))) std(abs_z_lesion_finalset(strcmp(Histo, TYPE))) ];
                    TYPE='IIA'; FCDTypeIIA_z_lesion_finalset(count, :) = [ mean(abs_z_lesion_finalset(strcmp(Histo, TYPE))) std(abs_z_lesion_finalset(strcmp(Histo, TYPE))) ];
                    
                    abs_z_flesion_finalset = abs(z_flesion_finalset(:, count));
                    TYPE='IIB'; FCDTypeIIB_z_flesion_finalset(count, :) = [ mean(abs_z_flesion_finalset(strcmp(Histo, TYPE))) std(abs_z_flesion_finalset(strcmp(Histo, TYPE))) ];
                    TYPE='IIA'; FCDTypeIIA_z_flesion_finalset(count, :) = [ mean(abs_z_flesion_finalset(strcmp(Histo, TYPE))) std(abs_z_flesion_finalset(strcmp(Histo, TYPE))) ];
                    count = count + 1;
                end
                
                mean_z_lesion_finalset = [];
                for i = 1 : num_of_FCD
                    mean_z_lesion_finalset(i, :) = mean(z_lesion_finalset(i, ~isnan(z_lesion_finalset(i, :))), 2);
                end
                
                p_type_IIa = [];
                p_type_IIb = [];
                p_type_IIa_vs_IIb = [];
                abs_z_flesion_finalset = flip(abs(z_flesion_finalset), 2);
                abs_z_fcont_finalset = flip(abs(squeeze(mean(z_fcont_finalset, 1))), 2);
                for i = 1 : size(abs_z_flesion_finalset, 2)
                    
                    mean_z_lesion_temp_a = mean(abs_z_flesion_finalset(strcmp(Histo, 'IIA'), i), 1);
                    std_z_lesion_temp_a  = std(abs_z_flesion_finalset(strcmp(Histo, 'IIA'), i), 0, 1);
                    mean_z_lesion_temp_cont = mean(abs_z_fcont_finalset(:, i), 1);
                    std_z_lesion_temp_cont  = std(abs_z_fcont_finalset(:, i), 0, 1);
                    pooled_std = sqrt(power(std_z_lesion_temp_a, 2)/num_of_FCDIIA+power(std_z_lesion_temp_cont, 2)/num_of_cont);
                    t = (mean_z_lesion_temp_a - mean_z_lesion_temp_cont)/pooled_std;
                    df = power(power(std_z_lesion_temp_a, 2)/num_of_FCDIIA+power(mean_z_lesion_temp_cont, 2)/num_of_cont, 2)/...
                        (power(power(std_z_lesion_temp_a, 2)/num_of_FCDIIA,2)/(num_of_FCDIIA-1) + power(power(mean_z_lesion_temp_cont, 2)/num_of_cont,2)/(num_of_cont-1));
                    p = tcdf(t, df);
                    p_type_IIa = [ p_type_IIa 2*min(p, 1-p) ];
                    
                    mean_z_lesion_temp_b = mean(abs_z_flesion_finalset(strcmp(Histo, 'IIB'), i), 1);
                    std_z_lesion_temp_b  = std(abs_z_flesion_finalset(strcmp(Histo, 'IIB'), i), 0, 1);
                    mean_z_lesion_temp_cont = mean(abs_z_fcont_finalset(:, i), 1);
                    std_z_lesion_temp_cont  = std(abs_z_fcont_finalset(:, i), 0, 1);
                    pooled_std = sqrt(power(std_z_lesion_temp_b, 2)/num_of_FCDIIB+power(std_z_lesion_temp_cont, 2)/num_of_cont);
                    t = (mean_z_lesion_temp_b - mean_z_lesion_temp_cont)/pooled_std;
                    df = power(power(std_z_lesion_temp_b, 2)/num_of_FCDIIB+power(mean_z_lesion_temp_cont, 2)/num_of_cont, 2)/...
                        (power(power(std_z_lesion_temp_b, 2)/num_of_FCDIIB,2)/(num_of_FCDIIB-1) + power(power(mean_z_lesion_temp_cont, 2)/num_of_cont,2)/(num_of_cont-1));
                    p = tcdf(t, df);
                    p_type_IIb = [ p_type_IIb 2*min(p, 1-p) ];
                    
                    pooled_std = sqrt(power(std_z_lesion_temp_a, 2)/num_of_FCDIIA+power(std_z_lesion_temp_b, 2)/num_of_FCDIIB);
                    t = (mean_z_lesion_temp_a - mean_z_lesion_temp_b)/pooled_std;
                    df = power(power(std_z_lesion_temp_a, 2)/num_of_FCDIIA+power(std_z_lesion_temp_b, 2)/num_of_FCDIIB, 2)/...
                        (power(power(std_z_lesion_temp_a, 2)/num_of_FCDIIA,2)/(num_of_FCDIIA-1) + power(power(std_z_lesion_temp_b, 2)/num_of_FCDIIB,2)/(num_of_FCDIIB-1));
                    p = tcdf(t, df);
                    p_type_IIa_vs_IIb = [ p_type_IIa_vs_IIb 2*min(p, 1-p) ];
                    
                end
                
                FDR_threshold = FDR([ p_type_IIa(1:10:51) p_type_IIb(1:10:51) p_type_IIa_vs_IIb(1:10:51) ], 0.05);
            end
            
            for visualize_result = 1
                abs_z_flesion_finalset = flip(abs(z_flesion_finalset), 2);
                linecolor = [0.8 0.8 0.8]; xtick = 1:10:51;
                
                f = figure; hold on;
                IIA_mem = find(strcmp(Histo, 'IIA')); IIA_mem(8) = [];
                IIB_mem = find(strcmp(Histo, 'IIB')); IIB_mem(8) = [];
                h1 = notBoxPlot(abs_z_flesion_finalset(strcmp(Histo, 'IIA'), xtick), xtick-1.5, 0.5, 'line');
                h2 = notBoxPlot(abs_z_flesion_finalset(strcmp(Histo, 'IIB'), xtick), xtick+1.5, 0.5, 'line');
                
                flag_mu  = 'off';
                flag_sd  = 'off';
                flag_sem = 'off';
                barhead_linewidth = 0.6;
                for i = 1 : size(h1, 2)
                    set(h1(i).data, 'markerfacecolor',[1,1,1],'color',[1,1,1]);
                    set(h2(i).data, 'markerfacecolor',[1,1,1],'color',[1,1,1]);
                    set(h1(i).mu,   'visible' , flag_mu, 'markerfacecolor',[1,1,1],'color',[1,1,1]);
                    set(h1(i).sd,   'visible' , flag_sd, 'LineStyle', ':', 'Color',[1,1,1]);
                    set(h1(i).sem,  'visible' , flag_sem, 'LineStyle', '-', 'Color',[1,1,1]);
                    set(h2(i).mu,   'visible' , flag_mu, 'markerfacecolor',[1,1,1],'color',[1,1,1]);
                    set(h2(i).sd,   'visible' , flag_sd, 'LineStyle', ':', 'Color',[1,1,1]);
                    set(h2(i).sem,  'visible' , flag_sem, 'LineStyle', '-', 'Color',[1,1,1]);
                    
                    x_1_temp = get(h1(i).mu, 'xData');
                    x_2_temp = get(h2(i).mu, 'xData');
                    mu_1 = get(h1(i).mu, 'YData');
                    mu_2 = get(h2(i).mu, 'YData');
                    sd_1 = get(h1(i).sd, 'YData');
                    sd_2 = get(h2(i).sd, 'YData');
                    med_1 = median(get(h1(i).data, 'YData'));
                    med_2 = median(get(h2(i).data, 'YData'));
                    
                    hold on; line([x_1_temp-barhead_linewidth x_1_temp+barhead_linewidth], [sd_1(1) sd_1(1)], 'LineWidth',1.3,'Color',[1 0 0]);
                    hold on; line([x_1_temp-barhead_linewidth x_1_temp+barhead_linewidth], [sd_1(2) sd_1(2)], 'LineWidth',1.3,'Color',[1 0 0]);
                    hold on; line([x_2_temp-barhead_linewidth x_2_temp+barhead_linewidth], [sd_2(1) sd_2(1)], 'LineWidth',1.3,'Color',[0 0 0]);
                    hold on; line([x_2_temp-barhead_linewidth x_2_temp+barhead_linewidth], [sd_2(2) sd_2(2)], 'LineWidth',1.3,'Color',[0 0 0]);
                    hold on; line([x_1_temp-barhead_linewidth x_1_temp+barhead_linewidth], [mu_1 mu_1], 'LineWidth',1.5,'Color',[1 0 0]);
                    hold on; line([x_2_temp-barhead_linewidth x_2_temp+barhead_linewidth], [mu_2 mu_2], 'LineWidth',1.5,'Color',[0 0 0]);
                    hold on; scatter(x_1_temp, mu_1, 30, 'MarkerEdgeColor', [1 0 0], 'markerfacecolor', [1 0 0]);
                    hold on; scatter(x_2_temp, mu_2, 30, 'MarkerEdgeColor', [0 0 0], 'markerfacecolor', [0 0 0]);
                end
                
                for i = 1 : length(IIA_mem)
                    data = smooth(abs_z_flesion_finalset(IIA_mem(i), :), 3);
                    plot((1:size(abs_z_flesion_finalset, 2))-1.5, data, 'Color', [1 0.8 0.8]);
                end
                for i = 1 : length(IIB_mem)
                    data = smooth(abs_z_flesion_finalset(IIB_mem(i), :), 3);
                    plot((1:size(abs_z_flesion_finalset, 2))+1.5, data, 'Color', [0.8 0.8 0.8]);
                end                
                
                plot(((1:10:length(r_thres_set))+1.5)', flip(FCDTypeIIB_z_flesion_finalset(1:10:51, 1), 1), 'k', 'LineWidth',1.5);
                plot(((1:10:length(r_thres_set))-1.5)', flip(FCDTypeIIA_z_flesion_finalset(1:10:51, 1), 1), 'r', 'LineWidth',1.5);
                
                set(gcf, 'Position', [587 500 600 400]);
                hp1 = plot([-5 70], [0 0]);
                set(hp1, 'Color', linecolor);
                set(gcf, 'Color', 'w');
                
                
                for i = 1 : size(h1, 2)
                    set(h1(i).mu,   'visible' , flag_mu, 'markerfacecolor',[1,1,1],'color',[1,1,1]);
                    set(h1(i).sd,   'visible' , flag_sd, 'LineStyle', ':', 'Color',[1,1,1]);
                    set(h1(i).sem,  'visible' , flag_sem, 'LineStyle', '-', 'Color',[1,1,1]);
                    set(h2(i).mu,   'visible' , flag_mu, 'markerfacecolor',[1,1,1],'color',[1,1,1]);
                    set(h2(i).sd,   'visible' , flag_sd, 'LineStyle', ':', 'Color',[1,1,1]);
                    set(h2(i).sem,  'visible' , flag_sem, 'LineStyle', '-', 'Color',[1,1,1]);
                    
                    x_1_temp = get(h1(i).mu, 'xData');
                    x_2_temp = get(h2(i).mu, 'xData');
                    mu_1 = get(h1(i).mu, 'YData');
                    mu_2 = get(h2(i).mu, 'YData');
                    sd_1 = get(h1(i).sd, 'YData');
                    sd_2 = get(h2(i).sd, 'YData');
                    med_1 = median(get(h1(i).data, 'YData'));
                    med_2 = median(get(h2(i).data, 'YData'));
                    
                    hold on; line([x_1_temp-barhead_linewidth x_1_temp+barhead_linewidth], [sd_1(1) sd_1(1)], 'LineWidth',1.3,'Color',[1 0 0]);
                    hold on; line([x_1_temp-barhead_linewidth x_1_temp+barhead_linewidth], [sd_1(2) sd_1(2)], 'LineWidth',1.3,'Color',[1 0 0]);
                    hold on; line([x_2_temp-barhead_linewidth x_2_temp+barhead_linewidth], [sd_2(1) sd_2(1)], 'LineWidth',1.3,'Color',[0 0 0]);
                    hold on; line([x_2_temp-barhead_linewidth x_2_temp+barhead_linewidth], [sd_2(2) sd_2(2)], 'LineWidth',1.3,'Color',[0 0 0]);
                    hold on; line([x_1_temp-barhead_linewidth x_1_temp+barhead_linewidth], [mu_1 mu_1], 'LineWidth',1.5,'Color',[1 0 0]);
                    hold on; line([x_2_temp-barhead_linewidth x_2_temp+barhead_linewidth], [mu_2 mu_2], 'LineWidth',1.5,'Color',[0 0 0]);
                    hold on; line([x_1_temp x_1_temp], [sd_1(1) sd_1(2)], 'LineWidth',1.3,'Color',[1 0 0]);
                    hold on; line([x_2_temp x_2_temp], [sd_2(1) sd_2(2)], 'LineWidth',1.3,'Color',[0 0 0]);
                    hold on; scatter(x_1_temp, mu_1, 30, 'MarkerEdgeColor', [1 0 0], 'markerfacecolor', [1 0 0]);
                    hold on; scatter(x_2_temp, mu_2, 30, 'MarkerEdgeColor', [0 0 0], 'markerfacecolor', [0 0 0]);
                end
                
                xlim([ -5 55 ]); ylim([ -0.2 2.2 ]);
                
                for i = 1 : 10: size(abs_z_flesion_finalset, 2)
                    if(p_type_IIa(i) <= FDR_threshold)
                        scatter(i-1.5, -0.1, 'r*');
                    elseif(p_type_IIa(i) <= 0.05 & p_type_IIa(i) > FDR_threshold)
                        scatter(i-1.5, -0.1, 'ro');
                    end
                    
                    if(p_type_IIb(i) <= FDR_threshold)
                        scatter(i+1.5, -0.1, 'k*');
                    elseif(p_type_IIb(i) <= 0.05 & p_type_IIb(i) > FDR_threshold)
                        scatter(i+1.5, -0.1, 'ko');
                    end
                    
                end
                
                xlabel = num2cell([flip(r_thres_set(1:10:51), 2)]);
                set(gca, 'XTick', [1:10:51]);
                set(gca, 'XTickLabel', xlabel);
                
                ODIR = '/local_raid/seokjun/01_project/04_IntracorticalAnalysis/03_result/lesion_profile/';
                if(save_file)
                    export_fig([ODIR '/43_FCD_lesion_rs-functional_connectivity_profile_typeIIa_IIb_new' ], '-m4', '-png'); close(gcf);
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% 016 057 066 072 081 085 086 089 -> 066 072 were selected after all...
                i = find(strcmp(CODE, '072')); % 2
                    
                for i = 1 : num_of_FCD
                    
                    disp(['case ' CODE{i} ' start!']);
                    if(exist([ LDIR '/' PREFIX{i} '_' CODE{i} '_label_union_left_rsl.txt' ], 'file'))
                        lesion_label_data = SurfStatReadData( { ...
                            [ LDIR '/' PREFIX{i} '_' CODE{i} '_label_union_left_rsl.txt'], ...
                            [ LDIR '/' PREFIX{i} '_' CODE{i} '_label_union_right_rsl.txt' ] } );
                    end
                    lesion_label_data = lesion_label_data([1:downsample_vert 40963:40963+10241]);
                    
                    %% 1) average time series within the lesion at this given patient
                    %% 2) compute the connectivity of the lesion with the rest of the brain
                    %% 3) transform the connectivity to Fishers' z-score
                    lesion_ts = mean(ts_set(:, lesion_label_data ~= 0, i), 2);
                    lesion_conn = corr(lesion_ts, ts_set(:, :, i));
                    flesion_conn = 0.5*(log((1+lesion_conn)./(1-lesion_conn)));
                    
                    %% 1) average time series within the lesional area at each control
                    %% 2) compute the connectivity of the area corresponding to the lesion with the rest of the brain
                    %% 3) transform the connectivity to Fishers' z-score
                    cont_conn = single(zeros(num_of_cont, length(surf_lowres.coord)));
                    fcont_conn = single(zeros(num_of_cont, length(surf_lowres.coord)));
                    for j = 1 : num_of_cont
                        fprintf([ num2str(j) ' ']);
                        cont_lesion_ts = mean(ts_set(:, lesion_label_data ~= 0, num_of_FCD+j), 2);
                        cont_conn(j, :) = corr(cont_lesion_ts, ts_set(:, :, num_of_FCD+j));
                        fcont_conn(j, :) = 0.5*(log((1+cont_conn(j, :))./(1-cont_conn(j, :))));
                    end
                    fprintf('\n');
                    
                    %% save these connecitivity measures
                    lesion_conn_set{i} = lesion_conn;
                    cont_conn_set{i}   = cont_conn;
                    
                    %% 1) compute the mean and sd of the connectivity matrices across all controls
                    %% 2) compute the z-score of the the connectivity matrices in this patient
                    %% 3) 'average_z_lesion_conn' saves the 'z_lesion_conn' across iterations
                    %%    'average_mean_cont_conn' saves the 'mean_cont_conn' across iterations
                    %% 4) compute internal z-score within this case
                    %% ** we do above procedures for both r and Fishers' score
                    
                    %% r
                    mean_cont_conn = mean(cont_conn, 1);
                    std_cont_conn  = std(cont_conn, 0, 1);
                    z_lesion_conn = (lesion_conn - mean_cont_conn)./std_cont_conn;
                    
                    %% fishers' score
                    mean_fcont_conn = mean(fcont_conn, 1);
                    std_fcont_conn  = std(fcont_conn, 0, 1);
                    z_flesion_conn = (flesion_conn - mean_fcont_conn)./std_fcont_conn;
                    
                    sm_mean_fcont_conn = SurfStatSmooth(mean_cont_conn, surf_lowres, 2);   sm_mean_fcont_conn(isnan(sm_mean_fcont_conn)) = 0; sm_mean_fcont_conn(isinf(sm_mean_fcont_conn)) = 0;
                    sm_flesion_conn    = SurfStatSmooth(flesion_conn, surf_lowres, 2);     sm_flesion_conn(isnan(sm_flesion_conn)) = 0; sm_flesion_conn(isinf(sm_flesion_conn)) = 0;
                    sm_z_flesion_conn    = SurfStatSmooth(z_flesion_conn, surf_lowres, 2); sm_z_flesion_conn(isnan(sm_z_flesion_conn)) = 0; sm_z_flesion_conn(isinf(sm_z_flesion_conn)) = 0;
                    
                    ODIR = '/local_raid/seokjun/01_project/05_NoelRest/03_Result/01_DPARSF_processing/DPARSFA_scrub/02_without_mean_sig_regout/seed_based_analysis_lesion/';
                    figure; SurfStatView(lesion_label_data, surf_lowres); colormap(flip(blackblue));
                    export_fig([ODIR '/mcd_' CODE{i} '_lesion_label' ], '-m4', '-png'); close(gcf);
                    
                    figure; SurfStatView(mean_fcont_conn, surf_lowres); SurfStatColLim([-1 1]);
                    export_fig([ODIR '/mcd_' CODE{i} '_fcontrol_mean' ], '-m4', '-png'); close(gcf);
%                     
%                     figure; SurfStatView(mean_fcont_conn, surf_lowres); SurfStatColLim([0.4 1]); colormap(red);
%                     export_fig([ODIR '/mcd_' CODE{i} '_fcontrol_mean_thres' ], '-m4', '-png'); close(gcf);
%                     
%                     figure; SurfStatView((mean_fcont_conn<0.5)*2, surf_lowres); SurfStatColLim([0 1]); colormap(blackblue);
%                     export_fig([ODIR '/mcd_' CODE{i} '_fcontrol_mean_binarized' ], '-m4', '-png'); close(gcf);
                    
                    figure; SurfStatView(sm_flesion_conn, surf_lowres); SurfStatColLim([-1 1]);
                    export_fig([ODIR '/mcd_' CODE{i} '_flesion' ], '-m4', '-png'); close(gcf);
                    
                    figure; SurfStatView(sm_z_flesion_conn, surf_lowres); SurfStatColLim([-3 3]);
                    export_fig([ODIR '/mcd_' CODE{i} '_zscore_flesion' ], '-m4', '-png'); close(gcf);
                    
                end
                
            end
           
        end
        for method_3 = 1
            
            %% data read and compute the connectivity
            lesion_conn_set                 = cell(num_of_FCD, 1);
            cont_conn_set                   = cell(num_of_FCD, 1);                        
            
            average_z_lesion_conn           = zeros(num_of_FCD, downsample_vert*2);
            internal_z_lesion_conn          = zeros(num_of_FCD, downsample_vert*2);
            average_mean_cont_conn          = zeros(num_of_FCD, downsample_vert*2);            
            
            average_z_flesion_conn          = zeros(num_of_FCD, downsample_vert*2);
            internal_z_flesion_conn         = zeros(num_of_FCD, downsample_vert*2);
            average_mean_fcont_conn         = zeros(num_of_FCD, downsample_vert*2);
            
            z_cont_conn_set                 = zeros(num_of_cont, downsample_vert*2, num_of_FCD);
            z_fcont_conn_set                = zeros(num_of_cont, downsample_vert*2, num_of_FCD);
            lesion_label_data_set           = zeros(num_of_FCD, downsample_vert*2);
            
            asso_conn_set                   = cell(num_of_FCD+num_of_cont, 1);
            asso_conn_pset                  = cell(num_of_FCD+num_of_cont, 1);
            
            parfor i = 1 : num_of_FCD+num_of_cont
                
                disp(['case ' CODE{i} ' start!']);
                if(exist([ LDIR '/' PREFIX{i} '_' CODE{i} '_label_union_left_rsl.txt' ], 'file'))
                    lesion_label_data = SurfStatReadData( { ...
                        [ LDIR '/' PREFIX{i} '_' CODE{i} '_label_union_left_rsl.txt'], ...
                        [ LDIR '/' PREFIX{i} '_' CODE{i} '_label_union_right_rsl.txt' ] } );
                    
                    lesion_label_data = lesion_label_data([1:downsample_vert 40963:40963+downsample_vert-1]);
                    lesion_label_data_set(i, :) = lesion_label_data;
                end
                %% 1) average time series within the lesion at this given patient
                %% 2) compute the connectivity of the lesion with the rest of the brain
                %% 3) transform the connectivity to Fishers' z-score                
                [ asso_conn asso_conn_p ] = corr(ts_set(:, :, i));
                asso_conn_set(i) = { single(asso_conn) };
                asso_conn_pset(i)= { single(asso_conn_p) };
                
            end
            
            asso_conn_set_temp = zeros(downsample_vert*2, downsample_vert*2, num_of_FCD+num_of_cont);
            asso_conn_pset_temp = zeros(downsample_vert*2, downsample_vert*2, num_of_FCD+num_of_cont);
            for i = 1 : num_of_FCD+num_of_cont
                asso_conn_set_temp(:, :, i) = asso_conn_set{i};
                asso_conn_pset_temp(:, :, i) = asso_conn_pset{i};
            end
            
            asso_conn_set = single(asso_conn_set_temp);
            asso_conn_pset = single(asso_conn_pset_temp);
            clear('asso_conn_set_temp', 'asso_conn_pset_temp');
            
            q_val = 0.001;
            thres = 0.05;
            n_step = 9;
            stepwise_conn_cont_set = zeros(num_of_cont, downsample_vert*2, n_step, num_of_FCD);
            normalized_stepwise_conn_cont_set = zeros(num_of_cont, downsample_vert*2, n_step, num_of_FCD);
            stepwise_conn_cont_p_set = zeros(n_step, downsample_vert*2, num_of_FCD);
            for i = 1 : num_of_FCD
                
                stepwise_conn_cont = cell(num_of_cont, 1);
                lesion_label_data = lesion_label_data_set(i, :);
                parfor j = 1 : num_of_cont
                    
                    j
                    asso_conn   = asso_conn_set(:, :, num_of_FCD+j);
                    asso_conn_p = asso_conn_pset(:, :, num_of_FCD+j); 
                    basso_conn_FDR = gretna_R2b(asso_conn.*(asso_conn>0), 's', thres);
                    
                    Wq = findwalks_npath(basso_conn_FDR, n_step); 
                    
                    temp = squeeze(sum(Wq(lesion_label_data>0, :, :), 1));
                    stepwise_conn_cont(j) = { temp' };
                    
                end
                
                temp = stepwise_conn_cont;
                stepwise_conn_cont = zeros(num_of_cont, downsample_vert*2, n_step);
                for j = 1 : n_step 
                    
                    for c = 1 : num_of_cont
                    
                        stepwise_conn_cont(c, :, j) = temp{c}(j, :);
                        
                    end
                    
                end
                clear('temp');
                
%                 figure; SurfStatView(mean(stepwise_conn_cont(:, :, 7), 1).*mask_lowres, surf_lowres);
               
                normalized_stepwise_conn_cont = stepwise_conn_cont;
                for s = 1 : n_step
                    for c = 1 : num_of_cont
                        
                        normalized_stepwise_conn_cont(c, :, s) = stepwise_conn_cont(c, :, s)/mean(stepwise_conn_cont(c, :, s));
                        
                    end
                end
                    
                stepwise_conn_cont_p = cell(n_step, 1);
                for s = 1 : n_step
                    stepwise_conn_cont_p_temp = cell(1, downsample_vert*2);
                    parfor j = 1 : downsample_vert*2                        
                        j                        
                        
                        [h, p, ci, stats] = ttest(normalized_stepwise_conn_cont(:, j, s));                        
                        stepwise_conn_cont_p_temp(j) = { p };
                    end
                    stepwise_conn_cont_p_temp = cell2mat(stepwise_conn_cont_p_temp);
                    stepwise_conn_cont_p_temp(isnan(stepwise_conn_cont_p_temp))=1;
                    stepwise_conn_cont_p(s) = { stepwise_conn_cont_p_temp };
                end
                
                temp = stepwise_conn_cont_p;
                stepwise_conn_cont_p = zeros(n_step, downsample_vert*2);
                for j = 1 : n_step 
                    
                    stepwise_conn_cont_p(j, :) = temp{j};                    
                    
                end
                clear('temp');                
                
%                 figure; SurfStatView(stepwise_conn_cont_p{5}, surf_lowres); SurfStatColLim([0 0.0001]); colormap(blue);
                
                stepwise_conn_cont_set(:, :, :, i) = stepwise_conn_cont;
                normalized_stepwise_conn_cont_set(:, :, :, i) = normalized_stepwise_conn_cont;
                stepwise_conn_cont_p_set(:, :, i) = stepwise_conn_cont_p;
            end     
        end        
    end
end

