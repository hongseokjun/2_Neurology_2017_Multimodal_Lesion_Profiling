clear all
close all
addpath('/data/noel/noel2/local/matlab/surfstat_chicago');

%% configure analysis parameters
meansig_cov = 'wbsignocov'; % wbsigcov or wbsignocov
smooth_method = 'minc';   % minc or surfstat
printfigs = 0;            % 1 (print) or 0 (no print)
load_mat = 1;             % 1 (load pre-saved files) or 0 (newly read every files)

%% setup directories
for setup_directories = 1
    if(strcmp(meansig_cov, 'wbsigcov'))
        WDIR = ['/host/gypsy/local_raid/seokjun/01_project/05_NoelRest/03_Result/01_DPARSF_processing/DPARSFA_scrub/01_with_mean_sig_regout/'];
    else
        WDIR = ['/host/gypsy/local_raid/seokjun/01_project/05_NoelRest/03_Result/01_DPARSF_processing/DPARSFA_scrub/02_without_mean_sig_regout/'];
    end
    CDIR = '/local_raid/seokjun/01_project/05_NoelRest/00_CaseList/';
    SDIR = '/local_raid/seokjun/01_project/05_NoelRest/01_Analysis/';    
    RDIR = '/local_raid/seokjun/01_project/05_NoelRest/03_Result/02_functional_network_local/';
    RDIR_topo = '//local_raid/seokjun/01_project/05_NoelRest/03_Result/03_functional_network_topology/';
    OUTPATH = '/local_raid/seokjun/01_project/04_IntracorticalAnalysis/03_result/lesion_profile/';
end

%% read surfaces
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
        
   rois = SurfStatReadData1([ SDIR '/parcellation/surf_reg_model_left_parcellation_kmeans.txt' ]);
   rois = [rois, rois+584];  
   
   if(printfigs)
       f     = figure; BoSurfStatViewData(rois,isurf,'');
       colormap([jet;jet])
   end
   
end

%% setup colormap
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

%% read csv 
for readcsv = 1
    
    fid   = fopen([ CDIR 'data_matlab_func_org.csv' ]) ;
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

%% Read pre-calculated z-score of T1 and FLAIR in the lesion
for read_matfiles = 1
        
    load('p_value_set_lesion_profile_all_modalities.mat');
    load('p_value_set_lesion_profile_T1_FLAIR.mat');
    load('p_value_set_lesion_profile_DTI.mat');
    load('p_value_set_lesion_cortical_depth_profile_T1_FLAIR.mat');
    load('p_value_set_lesion_cortical_depth_profile_DTI.mat');
    load('/local_raid/seokjun/01_project/04_IntracorticalAnalysis/01_analysis/pipeline/postprocessing_zscore_visualization/colormap_noel1.mat');
    
end

%% read local data
for read_local_data = 1
    if load_mat == 0
        
        % load local data
        for loaddata5 = 1
            namesleft  = [];
            namesright = [];
            
            % ALFF
            ALFF            = zeros(length(GROUP),81924);
            ALFF5_minc      = zeros(length(GROUP),81924);
            ALFF20_minc     = zeros(length(GROUP),81924);
            
            for i = 1:length(GROUP)
                namesleft{i,1}   = [FOLDER{i} '/bbr_rest/' PREFIX{i} '_' CODE{i} '/' PREFIX{i} '_' CODE{i} '_zALFFMap_mid_left_81920'];
                namesright{i,1}  = [FOLDER{i} '/bbr_rest/' PREFIX{i} '_' CODE{i} '/' PREFIX{i} '_' CODE{i} '_zALFFMap_mid_right_81920'];
                try
                    ALFF(i,:)           = SurfStatReadData({[ namesleft{i} '_rsl.txt' ] , [ namesright{i} '_rsl.txt' ]} );
                    ALFF5_minc(i,:) 	= SurfStatReadData({[ namesleft{i} '_sm5_rsl.txt' ] , [ namesright{i} '_sm5_rsl.txt' ]} );
                    ALFF20_minc(i,:)    = SurfStatReadData({[ namesleft{i} '_sm20_rsl.txt' ] , [ namesright{i} '_sm20_rsl.txt' ]} );
                catch
                    disp([ PREFIX{i} '_' CODE{i} ]);
                end
            end
            
            ALFF5_surfstat  = SurfStatSmooth(ALFF,surf,5);
            ALFF20_surfstat = SurfStatSmooth(ALFF,surf,20);
            
            if(printfigs)
                f = figure;
                BoSurfStatViewData(mean(ALFF5_minc,1).*mask,isurf,'');
                BoSurfStatColLim([-0.6 0.6]);
                
                f = figure;
                BoSurfStatViewData(mean(ALFF20_minc,1).*mask,isurf,'');
                BoSurfStatColLim([-0.6 0.6]);
                
                f = figure;
                BoSurfStatViewData(mean(ALFF5_surfstat,1).*mask,isurf,'');
                BoSurfStatColLim([-0.6 0.6]);
                
                f = figure;
                BoSurfStatViewData(mean(ALFF20_surfstat,1).*mask,isurf,'');
                BoSurfStatColLim([-0.6 0.6]);
            end
            
            
            % fALFF
            fALFF            = zeros(length(GROUP),81924);
            fALFF5_minc      = zeros(length(GROUP),81924);
            fALFF20_minc     = zeros(length(GROUP),81924);
            
            for i = 1:length(GROUP)
                namesleft{i,1}   = [FOLDER{i} '/bbr_rest/' PREFIX{i} '_' CODE{i} '/' PREFIX{i} '_' CODE{i} '_zfALFFMap_mid_left_81920'];
                namesright{i,1}  = [FOLDER{i} '/bbr_rest/' PREFIX{i} '_' CODE{i} '/' PREFIX{i} '_' CODE{i} '_zfALFFMap_mid_right_81920'];
                try
                    fALFF(i,:)           = SurfStatReadData({[ namesleft{i} '_rsl.txt' ] , [ namesright{i} '_rsl.txt' ]} );
                    fALFF5_minc(i,:) 	 = SurfStatReadData({[ namesleft{i} '_sm5_rsl.txt' ] , [ namesright{i} '_sm5_rsl.txt' ]} );
                    fALFF20_minc(i,:)    = SurfStatReadData({[ namesleft{i} '_sm20_rsl.txt' ] , [ namesright{i} '_sm20_rsl.txt' ]} );
                catch
                    disp([ PREFIX{i} '_' CODE{i} ]);
                end
            end
            
            fALFF5_surfstat  = SurfStatSmooth(fALFF,surf,5);
            fALFF20_surfstat = SurfStatSmooth(fALFF,surf,20);
            
            if(printfigs)
                f = figure;
                BoSurfStatViewData(mean(fALFF5_minc,1).*mask,isurf,'');
                BoSurfStatColLim([-0.2 0.2]);
                
                f = figure;
                BoSurfStatViewData(mean(fALFF20_minc,1).*mask,isurf,'');
                BoSurfStatColLim([-0.2 0.2]);
                
                f = figure;
                BoSurfStatViewData(mean(fALFF5_surfstat,1).*mask,isurf,'');
                BoSurfStatColLim([-0.2 0.2]);
                
                f = figure;
                BoSurfStatViewData(mean(fALFF20_surfstat,1).*mask,isurf,'');
                BoSurfStatColLim([-0.2 0.2]);
            end
            
            % DC
            DC            = zeros(length(GROUP),81924);
            DC5_minc      = zeros(length(GROUP),81924);
            DC20_minc     = zeros(length(GROUP),81924);
            
            for i = 1:length(GROUP)
                namesleft{i,1}   = [FOLDER{i} '/bbr_rest/' PREFIX{i} '_' CODE{i} '/' PREFIX{i} '_' CODE{i} '_zDegreeCentrality_PositiveBinarizedSumBrainMap_mid_left_81920'];
                namesright{i,1}  = [FOLDER{i} '/bbr_rest/' PREFIX{i} '_' CODE{i} '/' PREFIX{i} '_' CODE{i} '_zDegreeCentrality_PositiveBinarizedSumBrainMap_mid_right_81920'];
                try
                    DC(i,:)           = SurfStatReadData({[ namesleft{i} '_rsl.txt' ] , [ namesright{i} '_rsl.txt' ]} );
                    DC5_minc(i,:)     = SurfStatReadData({[ namesleft{i} '_sm5_rsl.txt' ] , [ namesright{i} '_sm5_rsl.txt' ]} );
                    DC20_minc(i,:)    = SurfStatReadData({[ namesleft{i} '_sm20_rsl.txt' ] , [ namesright{i} '_sm20_rsl.txt' ]} );
                catch
                    disp([ PREFIX{i} '_' CODE{i} ]);
                end
            end
            
            DC5_surfstat  = SurfStatSmooth(DC,surf,5);
            DC20_surfstat = SurfStatSmooth(DC,surf,20);
            
            if(printfigs)
                f = figure;
                BoSurfStatViewData(mean(DC5_minc,1).*mask,isurf,'');
                BoSurfStatColLim([-0.5 0.5]);
                
                f = figure;
                BoSurfStatViewData(mean(DC20_minc,1).*mask,isurf,'');
                BoSurfStatColLim([-0.5 0.5]);
                
                f = figure;
                BoSurfStatViewData(mean(DC5_surfstat,1).*mask,isurf,'');
                BoSurfStatColLim([-0.5 0.5]);
                
                f = figure;
                BoSurfStatViewData(mean(DC20_surfstat,1).*mask,isurf,'');
                BoSurfStatColLim([-0.5 0.5]);
            end
            
            % DC2: PositiveWeighted
            DC2            = zeros(length(GROUP),81924);
            DC2_5_minc      = zeros(length(GROUP),81924);
            DC2_20_minc     = zeros(length(GROUP),81924);
            
            for i = 1:length(GROUP)
                namesleft{i,1}   = [FOLDER{i} '/bbr_rest/' PREFIX{i} '_' CODE{i} '/' PREFIX{i} '_' CODE{i} '_zDegreeCentrality_PositiveWeightedSumBrainMap_mid_left_81920'];
                namesright{i,1}  = [FOLDER{i} '/bbr_rest/' PREFIX{i} '_' CODE{i} '/' PREFIX{i} '_' CODE{i} '_zDegreeCentrality_PositiveWeightedSumBrainMap_mid_right_81920'];
                try
                    DC2(i,:)           = SurfStatReadData({[ namesleft{i} '_rsl.txt' ] , [ namesright{i} '_rsl.txt' ]} );
                    DC2_5_minc(i,:)     = SurfStatReadData({[ namesleft{i} '_sm5_rsl.txt' ] , [ namesright{i} '_sm5_rsl.txt' ]} );
                    DC2_20_minc(i,:)    = SurfStatReadData({[ namesleft{i} '_sm20_rsl.txt' ] , [ namesright{i} '_sm20_rsl.txt' ]} );
                catch
                    disp([ PREFIX{i} '_' CODE{i} ]);
                end
            end
            
            DC2_5_surfstat  = SurfStatSmooth(DC2,surf,5);
            DC2_20_surfstat = SurfStatSmooth(DC2,surf,20);
            
            if(printfigs)
                f = figure;
                BoSurfStatViewData(mean(DC2_5_minc,1).*mask,isurf,'');
                BoSurfStatColLim([-0.5 0.5]);
                
                f = figure;
                BoSurfStatViewData(mean(DC2_20_minc,1).*mask,isurf,'');
                BoSurfStatColLim([-0.5 0.5]);
                
                f = figure;
                BoSurfStatViewData(mean(DC2_5_surfstat,1).*mask,isurf,'');
                BoSurfStatColLim([-0.5 0.5]);
                
                f = figure;
                BoSurfStatViewData(mean(DC2_20_surfstat,1).*mask,isurf,'');
                BoSurfStatColLim([-0.5 0.5]);
            end  
            
            % DC3: Weighted (pos + neg)
            DC3            = zeros(length(GROUP),81924);
            DC3_5_minc      = zeros(length(GROUP),81924);
            DC3_20_minc     = zeros(length(GROUP),81924);
            
            for i = 1:length(GROUP)
                namesleft{i,1}   = [FOLDER{i} '/bbr_rest/' PREFIX{i} '_' CODE{i} '/' PREFIX{i} '_' CODE{i} '_zDegreeCentrality_WeightedSumBrainMap_mid_left_81920'];
                namesright{i,1}  = [FOLDER{i} '/bbr_rest/' PREFIX{i} '_' CODE{i} '/' PREFIX{i} '_' CODE{i} '_zDegreeCentrality_WeightedSumBrainMap_mid_right_81920'];
                try
                    DC3(i,:)           = SurfStatReadData({[ namesleft{i} '_rsl.txt' ] , [ namesright{i} '_rsl.txt' ]} );
                    DC3_5_minc(i,:)     = SurfStatReadData({[ namesleft{i} '_sm5_rsl.txt' ] , [ namesright{i} '_sm5_rsl.txt' ]} );
                    DC3_20_minc(i,:)    = SurfStatReadData({[ namesleft{i} '_sm20_rsl.txt' ] , [ namesright{i} '_sm20_rsl.txt' ]} );
                catch
                    disp([ PREFIX{i} '_' CODE{i} ]);
                end
            end
            
            DC3_5_surfstat  = SurfStatSmooth(DC3,surf,5);
            DC3_20_surfstat = SurfStatSmooth(DC3,surf,20);
            
            if(printfigs)
                f = figure;
                BoSurfStatViewData(mean(DC3_5_minc,1).*mask,isurf,'');
                BoSurfStatColLim([-0.5 0.5]);
                
                f = figure;
                BoSurfStatViewData(mean(DC3_20_minc,1).*mask,isurf,'');
                BoSurfStatColLim([-0.5 0.5]);
                
                f = figure;
                BoSurfStatViewData(mean(DC3_5_surfstat,1).*mask,isurf,'');
                BoSurfStatColLim([-0.5 0.5]);
                
                f = figure;
                BoSurfStatViewData(mean(DC3_20_surfstat,1).*mask,isurf,'');
                BoSurfStatColLim([-0.5 0.5]);
            end  
            
            % EC
            EC            = zeros(length(GROUP),81924);
            EC5_minc      = zeros(length(GROUP),81924);
            EC20_minc     = zeros(length(GROUP),81924);
            
            for i = 1:length(GROUP)
                namesleft{i,1}   = [FOLDER{i} '/bbr_rest/' PREFIX{i} '_' CODE{i} '/' PREFIX{i} '_' CODE{i} '_zEigenvectorCentralityMap_mid_left_81920'];
                namesright{i,1}  = [FOLDER{i} '/bbr_rest/' PREFIX{i} '_' CODE{i} '/' PREFIX{i} '_' CODE{i} '_zEigenvectorCentralityMap_mid_right_81920'];
                try
                    EC(i,:)           = SurfStatReadData({[ namesleft{i} '_rsl.txt' ] , [ namesright{i} '_rsl.txt' ]} );
                    EC5_minc(i,:)     = SurfStatReadData({[ namesleft{i} '_sm5_rsl.txt' ] , [ namesright{i} '_sm5_rsl.txt' ]} );
                    EC20_minc(i,:)    = SurfStatReadData({[ namesleft{i} '_sm20_rsl.txt' ] , [ namesright{i} '_sm20_rsl.txt' ]} );
                catch
                    disp([ PREFIX{i} '_' CODE{i} ]);
                end
            end
            
            EC5_surfstat  = SurfStatSmooth(EC,surf,5);
            EC20_surfstat = SurfStatSmooth(EC,surf,20);
            
            if(printfigs)
                f = figure;
                BoSurfStatViewData(mean(EC5_minc,1).*mask,isurf,'');
                BoSurfStatColLim([-0.5 0.5]);
                
                f = figure;
                BoSurfStatViewData(mean(EC20_minc,1).*mask,isurf,'');
                BoSurfStatColLim([-0.5 0.5]);
                
                f = figure;
                BoSurfStatViewData(mean(EC5_surfstat,1).*mask,isurf,'');
                BoSurfStatColLim([-0.5 0.5]);
                
                f = figure;
                BoSurfStatViewData(mean(EC20_surfstat,1).*mask,isurf,'');
                BoSurfStatColLim([-0.5 0.5]);
            end            
            
            % ReHo
            ReHo            = zeros(length(GROUP),81924);
            ReHo5_minc      = zeros(length(GROUP),81924);
            ReHo20_minc     = zeros(length(GROUP),81924);
            
            for i = 1:length(GROUP)
                namesleft{i,1}   = [FOLDER{i} '/bbr_rest/' PREFIX{i} '_' CODE{i} '/' PREFIX{i} '_' CODE{i} '_zReHoMap_mid_left_81920'];
                namesright{i,1}  = [FOLDER{i} '/bbr_rest/' PREFIX{i} '_' CODE{i} '/' PREFIX{i} '_' CODE{i} '_zReHoMap_mid_right_81920'];
                try
                    ReHo(i,:)           = SurfStatReadData({[ namesleft{i} '_rsl.txt' ] , [ namesright{i} '_rsl.txt' ]} );
                    ReHo5_minc(i,:)     = SurfStatReadData({[ namesleft{i} '_sm5_rsl.txt' ] , [ namesright{i} '_sm5_rsl.txt' ]} );
                    ReHo20_minc(i,:)    = SurfStatReadData({[ namesleft{i} '_sm20_rsl.txt' ] , [ namesright{i} '_sm20_rsl.txt' ]} );
                catch
                    disp([ PREFIX{i} '_' CODE{i} ]);
                end
            end
            
            ReHo5_surfstat  = SurfStatSmooth(ReHo,surf,5);
            ReHo20_surfstat = SurfStatSmooth(ReHo,surf,20);
            
            if(printfigs)
                f = figure;
                BoSurfStatViewData(mean(ReHo5_minc,1).*mask,isurf,'');
                BoSurfStatColLim([-2 2]);
                
                f = figure;
                BoSurfStatViewData(mean(ReHo20_minc,1).*mask,isurf,'');
                BoSurfStatColLim([-2 2]);
                
                f = figure;
                BoSurfStatViewData(mean(ReHo5_surfstat,1).*mask,isurf,'');
                BoSurfStatColLim([-2 2]);
                
                f = figure;
                BoSurfStatViewData(mean(ReHo20_surfstat,1).*mask,isurf,'');
                BoSurfStatColLim([-2 2]);
            end
            
            % VMHC
            VMHC            = zeros(length(GROUP),81924);
            VMHC5_minc      = zeros(length(GROUP),81924);
            VMHC20_minc     = zeros(length(GROUP),81924);
            
            for i = 1:length(GROUP)
                namesleft{i,1}   = [FOLDER{i} '/bbr_rest/' PREFIX{i} '_' CODE{i} '/' PREFIX{i} '_' CODE{i} '_zVMHCMap_mid_left_81920'];
                namesright{i,1}  = [FOLDER{i} '/bbr_rest/' PREFIX{i} '_' CODE{i} '/' PREFIX{i} '_' CODE{i} '_zVMHCMap_mid_right_81920'];
                try
                    VMHC(i,:)           = SurfStatReadData({[ namesleft{i} '_rsl.txt' ] , [ namesright{i} '_rsl.txt' ]} );
                    VMHC5_minc(i,:)     = SurfStatReadData({[ namesleft{i} '_sm5_rsl.txt' ] , [ namesright{i} '_sm5_rsl.txt' ]} );
                    VMHC20_minc(i,:)    = SurfStatReadData({[ namesleft{i} '_sm20_rsl.txt' ] , [ namesright{i} '_sm20_rsl.txt' ]} );
                catch
                    disp([ PREFIX{i} '_' CODE{i} ]);
                end
            end
            
            VMHC5_surfstat  = SurfStatSmooth(VMHC,surf,5);
            VMHC20_surfstat = SurfStatSmooth(VMHC,surf,20);
            
            if(printfigs)
                f = figure;
                BoSurfStatViewData(mean(VMHC5_minc,1).*mask,isurf,'');
                BoSurfStatColLim([-2 2]);
                
                f = figure;
                BoSurfStatViewData(mean(VMHC20_minc,1).*mask,isurf,'');
                BoSurfStatColLim([-2 2]);
                
                f = figure;
                BoSurfStatViewData(mean(VMHC5_surfstat,1).*mask,isurf,'');
                BoSurfStatColLim([-2 2]);
                
                f = figure;
                BoSurfStatViewData(mean(VMHC20_surfstat,1).*mask,isurf,'');
                BoSurfStatColLim([-2 2]);
            end
            
            % Thickness
            CT            = zeros(length(GROUP),81924);
            CT5_minc      = zeros(length(GROUP),81924);
            CT20_minc     = zeros(length(GROUP),81924);
            
            for i = 1:length(GROUP)
                namesleft{i,1}   = [ '//media/fdb7fc29-462a-45ab-83e5-e928dede54f0/seokjun/01_project/05_NoelRest/03_Result/00_DPARSF_processing/DPARSFA_scrub/03_morphometrics/' PREFIX{i} '_' CODE{i} '/' PREFIX{i} '_' CODE{i} '_native_rms_rsl_tlink'];
                namesright{i,1}  = [ '//media/fdb7fc29-462a-45ab-83e5-e928dede54f0/seokjun/01_project/05_NoelRest/03_Result/00_DPARSF_processing/DPARSFA_scrub/03_morphometrics/' PREFIX{i} '_' CODE{i} '/' PREFIX{i} '_' CODE{i} '_native_rms_rsl_tlink'];
                try
                    CT(i,:)           = SurfStatReadData({[ namesleft{i} '_1mm_left.txt' ] , [ namesright{i} '_1mm_right.txt' ]} );
                    CT5_minc(i,:)     = SurfStatReadData({[ namesleft{i} '_5mm_left.txt' ] , [ namesright{i} '_5mm_right.txt' ]} );
                    CT20_minc(i,:)    = SurfStatReadData({[ namesleft{i} '_20mm_left.txt' ] , [ namesright{i} '_20mm_right.txt' ]} );
                catch
                    disp([ PREFIX{i} '_' CODE{i} ]);
                end
            end
            
            CT5_surfstat  = SurfStatSmooth(CT,surf,5);
            CT20_surfstat = SurfStatSmooth(CT,surf,20);
            
            if(printfigs)
                f = figure;
                BoSurfStatViewData(mean(CT5_minc,1).*mask,isurf,'');
                BoSurfStatColLim([1 5]);
                
                f = figure;
                BoSurfStatViewData(mean(CT20_minc,1).*mask,isurf,'');
                BoSurfStatColLim([1 5]);
                
                f = figure;
                BoSurfStatViewData(mean(CT5_surfstat,1).*mask,isurf,'');
                BoSurfStatColLim([1 5]);
                
                f = figure;
                BoSurfStatViewData(mean(CT20_surfstat,1).*mask,isurf,'');
                BoSurfStatColLim([1 5]);
            end
            
        end
        
        ALFF  = single(ALFF);
        fALFF = single(fALFF);
        DC    = single(DC);
        DC2   = single(DC2);
        DC3   = single(DC3);
        EC    = single(EC);
        ReHo  = single(ReHo);
        VMHC  = single(VMHC);
        CT    = single(CT);
        
        ALFF5_minc    = single(ALFF5_minc);
        fALFF5_minc   = single(fALFF5_minc);
        DC5_minc      = single(DC5_minc);
        DC2_5_minc    = single(DC2_5_minc);
        DC3_5_minc    = single(DC3_5_minc);
        EC5_minc      = single(EC5_minc);
        ReHo5_minc    = single(ReHo5_minc);
        VMHC5_minc    = single(VMHC5_minc);
        CT5_minc      = single(CT5_minc);
        
        ALFF20_minc       = single(ALFF20_minc);
        fALFF20_minc      = single(fALFF20_minc);
        DC20_minc         = single(DC20_minc);
        DC2_20_minc       = single(DC2_20_minc);
        DC3_20_minc       = single(DC3_20_minc);
        EC20_minc         = single(EC20_minc);
        ReHo20_minc       = single(ReHo20_minc);
        VMHC20_minc       = single(VMHC20_minc);
        CT20_minc         = single(CT20_minc);
        
        ALFF5_surfstat    = single(ALFF5_surfstat);
        fALFF5_surfstat   = single(fALFF5_surfstat);
        DC5_surfstat      = single(DC5_surfstat);
        DC2_5_surfstat    = single(DC2_5_surfstat);
        DC3_5_surfstat    = single(DC3_5_surfstat);
        EC5_surfstat      = single(EC5_surfstat);
        ReHo5_surfstat    = single(ReHo5_surfstat);
        VMHC5_surfstat    = single(VMHC5_surfstat);
        CT5_surfstat      = single(CT5_surfstat);
        
        ALFF20_surfstat    = single(ALFF20_surfstat);
        fALFF20_surfstat   = single(fALFF20_surfstat);
        DC20_surfstat      = single(DC20_surfstat);
        DC2_20_surfstat    = single(DC2_20_surfstat);
        DC3_20_surfstat    = single(DC3_20_surfstat);
        EC20_surfstat      = single(EC20_surfstat);
        ReHo20_surfstat    = single(ReHo20_surfstat);
        VMHC20_surfstat    = single(VMHC20_surfstat);
        CT20_surfstat      = single(CT20_surfstat);
        
        if(strcmp(smooth_method, 'minc'))
            
            ALFF5       = ALFF5_minc;
            fALFF5      = fALFF5_minc;
            DC5         = DC5_minc;
            DC2_5       = DC2_5_minc;
            DC3_5       = DC3_5_minc;
            EC5         = EC5_minc;
            ReHo5       = ReHo5_minc;
            VMHC5       = VMHC5_minc;
            CT5         = CT5_minc;
            
            ALFF20      = ALFF20_minc;
            fALFF20     = fALFF20_minc;
            DC20        = DC20_minc;
            DC2_20      = DC2_20_minc;
            DC3_20      = DC3_20_minc;
            EC20        = EC20_minc;
            ReHo20      = ReHo20_minc;
            VMHC20      = VMHC20_minc;
            CT20        = CT20_minc;
            
        elseif(strcmp(smooth_method, 'surfstat'))
            
            ALFF5       = ALFF5_surfstat;
            fALFF5      = fALFF5_surfstat;
            DC5         = DC5_surfstat;
            DC2_5       = DC2_5_surfstat;
            DC3_5       = DC3_5_surfstat;
            EC5         = EC5_surfstat;
            ReHo5       = ReHo5_surfstat;
            VMHC5       = VMHC5_surfstat;
            CT5         = CT5_surfstat;
            
            ALFF20      = ALFF20_surfstat;
            fALFF20     = fALFF20_surfstat;
            DC20        = DC20_surfstat;
            DC2_20      = DC2_20_surfstat;
            DC3_20      = DC3_20_surfstat;
            EC20        = EC20_surfstat;
            ReHo20      = ReHo20_surfstat;
            VMHC20      = VMHC20_surfstat;
            CT20        = CT20_surfstat;
            
        end        
        
        save([RDIR '/01_localfeature_' smooth_method '_' meansig_cov '2.mat'], ...
            'ALFF',   'fALFF',   'DC',   'DC2',    'DC3',    'EC',   'ReHo',   'VMHC', ...
            'ALFF5',  'fALFF5',  'DC5',  'DC2_5',  'DC3_5',  'EC5',  'ReHo5',  'VMHC5', ...
            'ALFF20', 'fALFF20', 'DC20', 'DC2_20', 'DC3_20', 'EC20', 'ReHo20', 'VMHC20');
    else
        load([RDIR '/01_localfeature_' smooth_method '_' meansig_cov '2.mat']);
    end
end

%% caluclate zscores
for calculate_zscores = 1
    NC     = find(strcmp(GROUP,'NC'));
    n      = length(GROUP);
    mALFF5 = mean(ALFF5(NC,:),1);
    sALFF5 = std(ALFF5(NC,:),0,1);
    ZALFF5 = (ALFF5 - repmat(mALFF5,n,1)) ./ repmat(sALFF5,n,1);
    
    mfALFF5 = mean(fALFF5(NC,:),1);
    sfALFF5 = std(fALFF5(NC,:),0,1);
    ZfALFF5 = (fALFF5 - repmat(mfALFF5,n,1)) ./ repmat(sfALFF5,n,1);
    
    mReHo5 = mean(ReHo5(NC,:),1);
    sReHo5 = std(ReHo5(NC,:),0,1);
    ZReHo5 = (ReHo5 - repmat(mReHo5,n,1)) ./ repmat(sReHo5,n,1);
    
    mVMHC5 = mean(VMHC5(NC,:),1);
    sVMHC5 = std(VMHC5(NC,:),0,1);
    ZVMHC5 = (VMHC5 - repmat(mVMHC5,n,1)) ./ repmat(sVMHC5,n,1);
    
    mDC5   = mean(DC5(NC,:),1);
    sDC5   = std(DC5(NC,:),0,1);
    ZDC5   = (DC5 - repmat(mDC5,n,1)) ./ repmat(sDC5,n,1);
    
    mDC2_5   = mean(DC2_5(NC,:),1);
    sDC2_5   = std(DC2_5(NC,:),0,1);
    ZDC2_5   = (DC2_5 - repmat(mDC2_5,n,1)) ./ repmat(sDC2_5,n,1);    

    mDC3_5   = mean(DC3_5(NC,:),1);
    sDC3_5   = std(DC3_5(NC,:),0,1);
    ZDC3_5   = (DC3_5 - repmat(mDC3_5,n,1)) ./ repmat(sDC3_5,n,1);    
    
    mEC5   = mean(EC5(NC,:),1);
    sEC5   = std(EC5(NC,:),0,1);
    ZEC5   = (EC5 - repmat(mEC5,n,1)) ./ repmat(sEC5,n,1);    
    
    for i = 1 : length(GROUP(strcmp(GROUP, 'FCD')))
        if(printfigs)
            f = figure;BoSurfStatViewData(ZALFF5(i,:).*mask,isurf,['zALFF ' CODE{i}]);
            colormap(red)
            BoSurfStatColLim([2 4])
            exportfigbo(f,[RPATH 'lesdet_' CODE{i} '_zALFF_pos.png'],'png',10);
            close(f);
            
            f = figure;BoSurfStatViewData(ZALFF5(i,:).*mask,isurf,['zALFF ' CODE{i}]);
            colormap(blackblue)
            BoSurfStatColLim([-4 -2])
            exportfigbo(f,[RPATH 'lesdet_' CODE{i} '_zALFF_neg.png'],'png',10);
            close(f);
            
            f = figure;BoSurfStatViewData(ZfALFF5(i,:).*mask,isurf,['zfALFF ' CODE{i}]);
            colormap(red)
            BoSurfStatColLim([2 3])
            exportfigbo(f,[RPATH 'lesdet_' CODE{i} '_zfALFF_pos.png'],'png',10);
            close(f);
            
            f = figure;BoSurfStatViewData(ZfALFF5(i,:).*mask,isurf,['zfALFF ' CODE{i}]);
            colormap(blackblue)
            BoSurfStatColLim([-3 -2])
            exportfigbo(f,[RPATH 'lesdet_' CODE{i} '_zfALFF_neg.png'],'png',10);
            close(f);
            
            f = figure;BoSurfStatViewData(ZReHo5(i,:).*mask,isurf,['zReHo ' CODE{i}]);
            colormap(red)
            BoSurfStatColLim([1 3])
            exportfigbo(f,[RPATH 'lesdet_' CODE{i} '_zReHo5_pos.png'],'png',10);
            close(f);
            f = figure;BoSurfStatViewData(ZReHo5(i,:).*mask,isurf,['zReHo ' CODE{i}]);
            colormap(blackblue)
            BoSurfStatColLim([-3 -1])
            exportfigbo(f,[RPATH 'lesdet_' CODE{i} '_zReHo5_neg.png'],'png',10);
            
            f = figure;SurfStatViewData(ZVMHC5(i,:).*mask,isurf,['zVMHC ' CODE{i}]);
            colormap(red)
            BoSurfStatColLim([2 4])
            exportfigbo(f,[RPATH 'lesdet_' CODE{i} '_zVMHC_pos.png'],'png',10);
            close(f);
        end
    end
end

%% lesion profiling
for lesion_profiling = 1
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% lesion feature extraction
    Lesion_ALFF5  = zeros(num_of_FCD, 1);
    Lesion_fALFF5 = zeros(num_of_FCD, 1);
    Lesion_ReHo5  = zeros(num_of_FCD, 1);
    Lesion_DC5    = zeros(num_of_FCD, 1);
    Lesion_DC2_5  = zeros(num_of_FCD, 1);
    Lesion_DC3_5  = zeros(num_of_FCD, 1);
    Lesion_EC5    = zeros(num_of_FCD, 1);
    Lesion_VMHC5  = zeros(num_of_FCD, 1);
    
    Lesion_ALFF5_cont  = zeros(num_of_cont, num_of_FCD);
    Lesion_fALFF5_cont = zeros(num_of_cont, num_of_FCD);
    Lesion_ReHo5_cont  = zeros(num_of_cont, num_of_FCD);
    Lesion_DC5_cont    = zeros(num_of_cont, num_of_FCD);
    Lesion_DC2_5_cont  = zeros(num_of_cont, num_of_FCD);
    Lesion_DC3_5_cont  = zeros(num_of_cont, num_of_FCD);
    Lesion_EC5_cont    = zeros(num_of_cont, num_of_FCD);
    Lesion_VMHC5_cont  = zeros(num_of_cont, num_of_FCD);
    
    %% lesion labels
    for readlesionlabels = 1
        
        LESIONS = zeros(num_of_FCD,81924);
        for i = 1 : length(GROUP(strcmp(GROUP, 'FCD')))
            namesleft{i,1}   = ['/local_raid/seokjun/CTFCD-1.2.0_64/Lesion/lesion_surf/' PREFIX{i} '_' CODE{i} '_label_union_left_rsl.txt'];
            namesright{i,1}  = ['/local_raid/seokjun/CTFCD-1.2.0_64/Lesion/lesion_surf/' PREFIX{i} '_' CODE{i} '_label_union_right_rsl.txt'];
            try
                LESIONS(i,:)    = SurfStatReadData({namesleft{i},namesright{i}});
                les_vert = find(LESIONS(i,:)~=0);                
                
                % lesion
                z = ZALFF5(i,les_vert);  Lesion_ALFF5(i)  = mean(ZALFF5(i,les_vert));
                z = ZfALFF5(i,les_vert); Lesion_fALFF5(i) = mean(ZfALFF5(i,les_vert));
                z = ZReHo5(i,les_vert);  Lesion_ReHo5(i)  = mean(ZReHo5(i,les_vert));
                z = ZDC5(i,les_vert);    Lesion_DC5(i)    = mean(ZDC5(i,les_vert));
                z = ZDC2_5(i,les_vert);  Lesion_DC2_5(i)  = mean(ZDC2_5(i,les_vert));
                z = ZDC3_5(i,les_vert);  Lesion_DC3_5(i)  = mean(ZDC3_5(i,les_vert));
                z = ZEC5(i,les_vert);    Lesion_EC5(i)    = mean(ZEC5(i,les_vert));
                z = ZVMHC5(i,les_vert);  Lesion_VMHC5(i)  = mean(ZVMHC5(i,les_vert));
                
                % healthy cortex
                z_ALFF5    = ZALFF5(strcmp(GROUP, 'NC'),les_vert);
                z_fALFF5   = ZfALFF5(strcmp(GROUP, 'NC'),les_vert);
                z_ReHo5    = ZReHo5(strcmp(GROUP, 'NC'),les_vert);
                z_DC5      = ZDC5(strcmp(GROUP, 'NC'),les_vert);
                z_DC2_5    = ZDC2_5(strcmp(GROUP, 'NC'),les_vert);
                z_DC3_5    = ZDC3_5(strcmp(GROUP, 'NC'),les_vert);
                z_EC5      = ZEC5(strcmp(GROUP, 'NC'),les_vert);
                z_VMHC5    = ZVMHC5(strcmp(GROUP, 'NC'),les_vert);
                
                for h = 1 : num_of_cont
                    z_ALFF5_ind = z_ALFF5(h, :);   Lesion_ALFF5_cont(h, i)  = mean(z_ALFF5_ind, 2);
                    z_fALFF5_ind = z_fALFF5(h, :); Lesion_fALFF5_cont(h, i) = mean(z_fALFF5_ind, 2);
                    z_ReHo5_ind = z_ReHo5(h, :);   Lesion_ReHo5_cont(h, i)  = mean(z_ReHo5_ind, 2);
                    z_DC5_ind = z_DC5(h, :);       Lesion_DC5_cont(h, i)    = mean(z_DC5_ind, 2);
                    z_DC2_5_ind = z_DC2_5(h, :);   Lesion_DC2_5_cont(h, i)  = mean(z_DC2_5_ind, 2);
                    z_DC3_5_ind = z_DC3_5(h, :);   Lesion_DC3_5_cont(h, i)  = mean(z_DC3_5_ind, 2);
                    z_EC5_ind = z_EC5(h, :);       Lesion_EC5_cont(h, i)    = mean(z_EC5_ind, 2);
                    z_VMHC5_ind = z_VMHC5(h, :);   Lesion_VMHC5_cont(h, i)  = mean(z_VMHC5_ind, 2);
                end                                
            
            catch
                disp(namesleft{i})
            end
        end
        
    end    
     
    % Parcellation-based centrality lesion profiling
    load([ RDIR_topo '/06_individual_centrality_' meansig_cov '.mat' ]);
    
    %% lesion labels
    for readlesionlabels = 1
        
        LESIONS = zeros(sum(strcmp(GROUP, 'FCD')),81924);
        for i = 1 : length(GROUP(strcmp(GROUP, 'FCD')))
            namesleft{i,1}   = ['/local_raid/seokjun/CTFCD-1.2.0_64/Lesion/lesion_surf/' PREFIX{i} '_' CODE{i} '_label_union_left_rsl.txt'];
            namesright{i,1}  = ['/local_raid/seokjun/CTFCD-1.2.0_64/Lesion/lesion_surf/' PREFIX{i} '_' CODE{i} '_label_union_right_rsl.txt'];
            try
                LESIONS(i,:)    = SurfStatReadData({namesleft{i},namesright{i}});
                les_vert = find(LESIONS(i,:)~=0);                
                
                % lesion
                z = zDC_ind_mapping(i,les_vert); included_vert = find(z < (mean(z)+2*std(z)) & z > (mean(z)-2*std(z)));
                Lesion_DC(i) = mean(zDC_ind_mapping(i, les_vert(included_vert)));
                z = zBC_ind_mapping(i,les_vert); included_vert = find(z < (mean(z)+2*std(z)) & z > (mean(z)-2*std(z)));
                Lesion_BC(i) = mean(zBC_ind_mapping(i,les_vert(included_vert)));
                z = zEC_ind_mapping(i,les_vert); included_vert = find(z < (mean(z)+2*std(z)) & z > (mean(z)-2*std(z)));
                Lesion_EC(i) = mean(zEC_ind_mapping(i,les_vert(included_vert)));
                z = zPC_ind_mapping(i,les_vert); included_vert = find(z < (mean(z)+2*std(z)) & z > (mean(z)-2*std(z)));
                Lesion_PC(i) = mean(zPC_ind_mapping(i,les_vert(included_vert)));
                
                % healthy cortex
                z_DC = zDC_ind_mapping(strcmp(GROUP, 'NC'),les_vert);
                z_BC = zDC_ind_mapping(strcmp(GROUP, 'NC'),les_vert);
                z_EC = zDC_ind_mapping(strcmp(GROUP, 'NC'),les_vert);
                z_PC = zDC_ind_mapping(strcmp(GROUP, 'NC'),les_vert);
                
                for h = 1 : sum(strcmp(GROUP, 'NC'))
                    z_DC_ind = z_DC(h, :); included_vert_DC = find(z_DC_ind < (mean(z_DC_ind)+2*std(z_DC_ind)) & z_DC_ind > (mean(z_DC_ind)-2*std(z_DC_ind)));
                    z_BC_ind = z_BC(h, :); included_vert_BC = find(z_BC_ind < (mean(z_BC_ind)+2*std(z_BC_ind)) & z_BC_ind > (mean(z_BC_ind)-2*std(z_BC_ind)));
                    z_EC_ind = z_EC(h, :); included_vert_EC = find(z_EC_ind < (mean(z_EC_ind)+2*std(z_EC_ind)) & z_EC_ind > (mean(z_EC_ind)-2*std(z_EC_ind)));
                    z_PC_ind = z_PC(h, :); included_vert_PC = find(z_PC_ind < (mean(z_PC_ind)+2*std(z_PC_ind)) & z_PC_ind > (mean(z_PC_ind)-2*std(z_PC_ind)));
                    
                    Lesion_DC_cont(h, i) = mean(zDC_ind_mapping(sum(strcmp(GROUP, 'FCD'))+h,les_vert(included_vert_DC)), 2);
                    Lesion_BC_cont(h, i) = mean(zBC_ind_mapping(sum(strcmp(GROUP, 'FCD'))+h,les_vert(included_vert_BC)), 2);
                    Lesion_EC_cont(h, i) = mean(zEC_ind_mapping(sum(strcmp(GROUP, 'FCD'))+h,les_vert(included_vert_EC)), 2);
                    Lesion_PC_cont(h, i) = mean(zPC_ind_mapping(sum(strcmp(GROUP, 'FCD'))+h,les_vert(included_vert_PC)), 2);
                end
            
%                 % lesion
%                 Lesion_DC(i) = mean(zDC_ind_mapping(i,les_vert));
%                 Lesion_BC(i) = mean(zBC_ind_mapping(i,les_vert));
%                 Lesion_EC(i) = mean(zEC_ind_mapping(i,les_vert));
%                 Lesion_PC(i) = mean(zPC_ind_mapping(i,les_vert));
%                 
%                 % healthy cortex
%                 Lesion_DC_cont(:, i) = mean(zDC_ind_mapping(strcmp(GROUP, 'NC'),les_vert), 2);
%                 Lesion_BC_cont(:, i) = mean(zBC_ind_mapping(strcmp(GROUP, 'NC'),les_vert), 2);
%                 Lesion_EC_cont(:, i) = mean(zEC_ind_mapping(strcmp(GROUP, 'NC'),les_vert), 2);
%                 Lesion_PC_cont(:, i) = mean(zPC_ind_mapping(strcmp(GROUP, 'NC'),les_vert), 2);
              
            catch
                disp(namesleft{i})
            end
        end
    end
    
    save_mat = 0;
    if(save_mat == 1)
        mean_z_lesion      = [ Lesion_ALFF5'; Lesion_ReHo5'; Lesion_DC; Lesion_BC ];
        mean_z_lesion_cont = cat(1, reshape(Lesion_ALFF5_cont', 1, size(Lesion_ALFF5_cont', 1), size(Lesion_ALFF5_cont', 2)), ...
                                    reshape(Lesion_ReHo5_cont', 1, size(Lesion_ReHo5_cont', 1), size(Lesion_ReHo5_cont', 2)), ...
                                    reshape(Lesion_DC_cont',    1, size(Lesion_DC_cont', 1),    size(Lesion_DC_cont', 2)),    ...
                                    reshape(Lesion_BC_cont',    1, size(Lesion_BC_cont', 1),    size(Lesion_BC_cont', 2)));
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    %% statistical tests
    for method_localfeat = 1
        FCD = find(strcmp(GROUP, 'FCD'));
        NC  = find(strcmp(GROUP, 'NC'));
        FCD_IIa = find(strcmp(Histo, 'IIA'))'; % FCD_IIa(FCD_IIa==16 | FCD_IIa==20) = []; 
        FCD_IIb = find(strcmp(Histo, 'IIB'))';
        
        p_type     = [];
        p_type_IIa = [];
        p_type_IIb = [];
        p_type_IIa_vs_IIb = [];
        p_type_ss     = [];             % surfstat
        p_type_IIa_ss = [];             % surfstat
        p_type_IIb_ss = [];             % surfstat
        p_type_IIa_vs_IIb_ss = [];      % surfstat

        
        %% feature configureation: ALFF5, fALFF5, ReHo5, DC5
        for feature_ALFF5 = 1
            feature_temp = Lesion_ALFF5;
            
            mean_z_lesion_temp = mean(feature_temp);
            std_z_lesion_temp  = std(feature_temp);
            pooled_std = sqrt(power(std_z_lesion_temp, 2)/num_of_FCD+1/num_of_cont);
            t = (mean_z_lesion_temp - 0)/pooled_std;
            df = power(power(std_z_lesion_temp, 2)/num_of_FCD+1/num_of_cont, 2)/(power(power(std_z_lesion_temp, 2)/num_of_FCD,2)/(num_of_FCD-1) + power(1/num_of_cont, 2)/(num_of_cont-1));
            p = tcdf(t, df);
            p_type = [ p_type 2*min(p, 1-p) ]
            g = ones(length(feature_temp), 1); M = 1 + term(g);
            slm = SurfStatLinMod(feature_temp, M);
            slm = SurfStatT(slm, g); p = 1 - tcdf(slm.t,slm.df);
            p_type_ss = [ p_type_ss 2*min(p, 1-p) ];
            
            mean_z_lesion_temp = mean(feature_temp(FCD_IIa));
            std_z_lesion_temp  = std(feature_temp(FCD_IIa));
            pooled_std = sqrt(power(std_z_lesion_temp, 2)/num_of_FCDIIA+1/num_of_cont);
            t = (mean_z_lesion_temp - 0)/pooled_std;
            df = power(power(std_z_lesion_temp, 2)/num_of_FCDIIA+1/num_of_cont, 2)/(power(power(std_z_lesion_temp, 2)/num_of_FCDIIA,2)/(num_of_FCDIIA-1) + power(1/num_of_cont, 2)/(num_of_cont-1));
            p = tcdf(t, df);
            p_type_IIa = [ p_type_IIa 2*min(p, 1-p) ]
            g = ones(length(feature_temp(FCD_IIa)), 1); M = 1 + term(g);
            slm = SurfStatLinMod(feature_temp(FCD_IIa), M);
            slm = SurfStatT(slm, g); p = 1 - tcdf(slm.t,slm.df);
            p_type_IIa_ss = [ p_type_IIa_ss 2*min(p, 1-p) ];
            
            mean_z_lesion_temp = mean(feature_temp(FCD_IIb));
            std_z_lesion_temp  = std(feature_temp(FCD_IIb));
            pooled_std = sqrt(power(std_z_lesion_temp, 2)/num_of_FCDIIB+1/num_of_cont);
            t = (mean_z_lesion_temp - 0)/pooled_std;
            df = power(power(std_z_lesion_temp, 2)/num_of_FCDIIB+1/num_of_cont, 2)/(power(power(std_z_lesion_temp, 2)/num_of_FCDIIB,2)/(num_of_FCDIIB-1) + power(1/num_of_cont, 2)/(num_of_cont-1));
            p = tcdf(t, df);
            p_type_IIb = [ p_type_IIb 2*min(p, 1-p) ]
            g = ones(length(feature_temp(FCD_IIb)), 1); M = 1 + term(g);
            slm = SurfStatLinMod(feature_temp(FCD_IIb), M);
            slm = SurfStatT(slm, g); p = 1 - tcdf(slm.t,slm.df);
            p_type_IIb_ss = [ p_type_IIb_ss 2*min(p, 1-p) ];
            
            mean_z_lesion_temp_a = mean(feature_temp(FCD_IIa));
            std_z_lesion_temp_a  = std(feature_temp(FCD_IIa));
            mean_z_lesion_temp_b = mean(feature_temp(FCD_IIb));
            std_z_lesion_temp_b  = std(feature_temp(FCD_IIb));
            pooled_std = sqrt(power(std_z_lesion_temp_a, 2)/num_of_FCDIIA+power(std_z_lesion_temp_b, 2)/num_of_FCDIIB);
            t = (mean_z_lesion_temp_a - mean_z_lesion_temp_b)/pooled_std;
            df = power(power(std_z_lesion_temp_a, 2)/num_of_FCDIIA+power(std_z_lesion_temp_b, 2)/num_of_FCDIIB, 2)/...
                (power(power(std_z_lesion_temp_a, 2)/num_of_FCDIIA,2)/(num_of_FCDIIA-1) + power(power(std_z_lesion_temp_b, 2)/num_of_FCDIIB,2)/(num_of_FCDIIB-1));
            p = tcdf(t, df);
            p_type_IIa_vs_IIb = [ p_type_IIa_vs_IIb min(p, 1-p) ]
            g = [-1*ones(length(feature_temp(FCD_IIa)), 1); ones(length(feature_temp(FCD_IIb)), 1); ]; M = 1 + term(g);
            slm = SurfStatLinMod([feature_temp(FCD_IIa); feature_temp(FCD_IIb)], M);
            slm = SurfStatT(slm, g); p = 1 - tcdf(slm.t,slm.df);
            p_type_IIa_vs_IIb_ss = [ p_type_IIa_vs_IIb_ss 2*min(p, 1-p) ];      
            p_value_set_ALFF = [ p_type_IIa p_type_IIb ];
            
            %% Kolmogorov-sminorv test for normal distribution
            % IIA
            for i = 1 : size(feature_temp(FCD_IIa), 2)
                
                xCDF = sort(feature_temp(FCD_IIa)); yCDF = normcdf(xCDF, mean(xCDF), std(xCDF));
                [h p] = kstest(feature_temp(FCD_IIa), 'CDF', [xCDF yCDF])
                
            end
            
            % IIB
            for i = 1 : size(feature_temp(FCD_IIb), 2)
                
                xCDF = sort(feature_temp(FCD_IIb)); yCDF = normcdf(xCDF, mean(xCDF), std(xCDF));
                [h p] = kstest(feature_temp(FCD_IIb), 'CDF', [xCDF yCDF])
                
            end
            
            %% effect size of FCD Type IIA
            pwer = 0.8;
            alph = 0.025;
            
            mean_z_lesion_temp = mean(feature_temp(FCD_IIa));
            std_z_lesion_temp  = std(feature_temp(FCD_IIa));
            pooled_std = sqrt(((num_of_FCDIIA-1)*power(std_z_lesion_temp, 2)+(num_of_cont-1)*1)/(num_of_FCDIIA+num_of_cont-2));
            effect_size_IIA = (mean_z_lesion_temp - 0)/pooled_std;
%             samplesize = sampsizepwr('t2',[0 1],[mean_z_lesion_temp],pwer,[],'Alpha', alph)
            
        end        
        for feature_ReHo5 = 1
            feature_temp = Lesion_ReHo5;
            
            mean_z_lesion_temp = mean(feature_temp);
            std_z_lesion_temp  = std(feature_temp);
            pooled_std = sqrt(power(std_z_lesion_temp, 2)/num_of_FCD+1/num_of_cont);
            t = (mean_z_lesion_temp - 0)/pooled_std;
            df = power(power(std_z_lesion_temp, 2)/num_of_FCD+1/num_of_cont, 2)/(power(power(std_z_lesion_temp, 2)/num_of_FCD,2)/(num_of_FCD-1) + power(1/num_of_cont, 2)/(num_of_cont-1));
            p = tcdf(t, df);
            p_type = [ p_type 2*min(p, 1-p) ]
            g = ones(length(feature_temp), 1); M = 1 + term(g);
            slm = SurfStatLinMod(feature_temp, M);
            slm = SurfStatT(slm, g); p = 1 - tcdf(slm.t,slm.df);
            p_type_ss = [ p_type_ss 2*min(p, 1-p) ];
            
            mean_z_lesion_temp = mean(feature_temp(FCD_IIa));
            std_z_lesion_temp  = std(feature_temp(FCD_IIa));
            pooled_std = sqrt(power(std_z_lesion_temp, 2)/num_of_FCDIIA+1/num_of_cont);
            t = (mean_z_lesion_temp - 0)/pooled_std;
            df = power(power(std_z_lesion_temp, 2)/num_of_FCDIIA+1/num_of_cont, 2)/(power(power(std_z_lesion_temp, 2)/num_of_FCDIIA,2)/(num_of_FCDIIA-1) + power(1/num_of_cont, 2)/(num_of_cont-1));
            p = tcdf(t, df);
            p_type_IIa = [ p_type_IIa 2*min(p, 1-p) ]
            g = ones(length(feature_temp(FCD_IIa)), 1); M = 1 + term(g);
            slm = SurfStatLinMod(feature_temp(FCD_IIa), M);
            slm = SurfStatT(slm, g); p = 1 - tcdf(slm.t,slm.df);
            p_type_IIa_ss = [ p_type_IIa_ss 2*min(p, 1-p) ];
            
            mean_z_lesion_temp = mean(feature_temp(FCD_IIb));
            std_z_lesion_temp  = std(feature_temp(FCD_IIb));
            pooled_std = sqrt(power(std_z_lesion_temp, 2)/num_of_FCDIIB+1/num_of_cont);
            t = (mean_z_lesion_temp - 0)/pooled_std;
            df = power(power(std_z_lesion_temp, 2)/num_of_FCDIIB+1/num_of_cont, 2)/(power(power(std_z_lesion_temp, 2)/num_of_FCDIIB,2)/(num_of_FCDIIB-1) + power(1/num_of_cont, 2)/(num_of_cont-1));
            p = tcdf(t, df);
            p_type_IIb = [ p_type_IIb 2*min(p, 1-p) ]
            g = ones(length(feature_temp(FCD_IIb)), 1); M = 1 + term(g);
            slm = SurfStatLinMod(feature_temp(FCD_IIb), M);
            slm = SurfStatT(slm, g); p = 1 - tcdf(slm.t,slm.df);
            p_type_IIb_ss = [ p_type_IIb_ss 2*min(p, 1-p) ];
            
            mean_z_lesion_temp_a = mean(feature_temp(FCD_IIa));
            std_z_lesion_temp_a  = std(feature_temp(FCD_IIa));
            mean_z_lesion_temp_b = mean(feature_temp(FCD_IIb));
            std_z_lesion_temp_b  = std(feature_temp(FCD_IIb));
            pooled_std = sqrt(power(std_z_lesion_temp_a, 2)/num_of_FCDIIA+power(std_z_lesion_temp_b, 2)/num_of_FCDIIB);
            t = (mean_z_lesion_temp_a - mean_z_lesion_temp_b)/pooled_std;
            df = power(power(std_z_lesion_temp_a, 2)/num_of_FCDIIA+power(std_z_lesion_temp_b, 2)/num_of_FCDIIB, 2)/...
                (power(power(std_z_lesion_temp_a, 2)/num_of_FCDIIA,2)/(num_of_FCDIIA-1) + power(power(std_z_lesion_temp_b, 2)/num_of_FCDIIB,2)/(num_of_FCDIIB-1));
            p = tcdf(t, df);
            p_type_IIa_vs_IIb = [ p_type_IIa_vs_IIb min(p, 1-p) ]
            g = [-1*ones(length(feature_temp(FCD_IIa)), 1); ones(length(feature_temp(FCD_IIb)), 1); ]; M = 1 + term(g);
            slm = SurfStatLinMod([feature_temp(FCD_IIa); feature_temp(FCD_IIb)], M);
            slm = SurfStatT(slm, g); p = 1 - tcdf(slm.t,slm.df);
            p_type_IIa_vs_IIb_ss = [ p_type_IIa_vs_IIb_ss 2*min(p, 1-p) ];
            
            p_value_set_ALFF_ReHo = [ p_type_IIa p_type_IIb p_type_IIa_vs_IIb ];
            p_value_set_ReHo = [ p_type_IIa p_type_IIb ];
            
            %% Kolmogorov-sminorv test for normal distribution
            % IIA
            for i = 1 : size(feature_temp(FCD_IIa), 2)
                
                xCDF = sort(feature_temp(FCD_IIa)); yCDF = normcdf(xCDF, mean(xCDF), std(xCDF));
                [h p] = kstest(feature_temp(FCD_IIa), 'CDF', [xCDF yCDF])
                
            end
            
            % IIB
            for i = 1 : size(feature_temp(FCD_IIb), 2)
                
                xCDF = sort(feature_temp(FCD_IIb)); yCDF = normcdf(xCDF, mean(xCDF), std(xCDF));
                [h p] = kstest(feature_temp(FCD_IIb), 'CDF', [xCDF yCDF])
                
            end            
            
            %% effect size of FCD Type IIA
            pwer = 0.8;
            alph = 0.025;
            
            mean_z_lesion_temp = mean(feature_temp(FCD_IIa));
            std_z_lesion_temp  = std(feature_temp(FCD_IIa));
            pooled_std = sqrt(((num_of_FCDIIA-1)*power(std_z_lesion_temp, 2)+(num_of_cont-1)*1)/(num_of_FCDIIA+num_of_cont-2));
            effect_size_IIA = (mean_z_lesion_temp - 0)/pooled_std;
%             samplesize = sampsizepwr('t2',[0 1],[mean_z_lesion_temp],pwer,[],'Alpha', alph)
            
        end
        for feature_DC = 1
            
            feature_temp = Lesion_DC;
            
            mean_z_lesion_temp = mean(feature_temp);
            std_z_lesion_temp  = std(feature_temp);
            pooled_std = sqrt(power(std_z_lesion_temp, 2)/num_of_FCD+1/num_of_cont);
            t = (mean_z_lesion_temp - 0)/pooled_std;
            df = power(power(std_z_lesion_temp, 2)/num_of_FCD+1/num_of_cont, 2)/(power(power(std_z_lesion_temp, 2)/num_of_FCD,2)/(num_of_FCD-1) + power(1/num_of_cont, 2)/(num_of_cont-1));
            p = tcdf(t, df);
            p_type = [ p_type 2*min(p, 1-p) ]
%             g = ones(length(feature_temp), 1); M = 1 + term(g);
%             slm = SurfStatLinMod(feature_temp, M);
%             slm = SurfStatT(slm, g); p = 1 - tcdf(slm.t,slm.df);
%             p_type_ss = [ p_type_ss 2*min(p, 1-p) ];
            
            mean_z_lesion_temp = mean(feature_temp(FCD_IIa));
            std_z_lesion_temp  = std(feature_temp(FCD_IIa));
            pooled_std = sqrt(power(std_z_lesion_temp, 2)/num_of_FCDIIA+1/num_of_cont);
            t = (mean_z_lesion_temp - 0)/pooled_std;
            df = power(power(std_z_lesion_temp, 2)/num_of_FCDIIA+1/num_of_cont, 2)/(power(power(std_z_lesion_temp, 2)/num_of_FCDIIA,2)/(num_of_FCDIIA-1) + power(1/num_of_cont, 2)/(num_of_cont-1));
            p = tcdf(t, df);
            p_type_IIa = [ p_type_IIa 2*min(p, 1-p) ]
%             g = ones(length(feature_temp(FCD_IIa)), 1); M = 1 + term(g);
%             slm = SurfStatLinMod(feature_temp(FCD_IIa), M);
%             slm = SurfStatT(slm, g); p = 1 - tcdf(slm.t,slm.df);
%             p_type_IIa_ss = [ p_type_IIa_ss 2*min(p, 1-p) ];
            
            mean_z_lesion_temp = mean(feature_temp(FCD_IIb));
            std_z_lesion_temp  = std(feature_temp(FCD_IIb));
            pooled_std = sqrt(power(std_z_lesion_temp, 2)/num_of_FCDIIB+1/num_of_cont);
            t = (mean_z_lesion_temp - 0)/pooled_std;
            df = power(power(std_z_lesion_temp, 2)/num_of_FCDIIB+1/num_of_cont, 2)/(power(power(std_z_lesion_temp, 2)/num_of_FCDIIB,2)/(num_of_FCDIIB-1) + power(1/num_of_cont, 2)/(num_of_cont-1));
            p = tcdf(t, df);
            p_type_IIb = [ p_type_IIb 2*min(p, 1-p) ]
%             g = ones(length(feature_temp(FCD_IIb)), 1); M = 1 + term(g);
%             slm = SurfStatLinMod(feature_temp(FCD_IIb), M);
%             slm = SurfStatT(slm, g); p = 1 - tcdf(slm.t,slm.df);
%             p_type_IIb_ss = [ p_type_IIb_ss 2*min(p, 1-p) ];
            
            mean_z_lesion_temp_a = mean(feature_temp(FCD_IIa));
            std_z_lesion_temp_a  = std(feature_temp(FCD_IIa));
            mean_z_lesion_temp_b = mean(feature_temp(FCD_IIb));
            std_z_lesion_temp_b  = std(feature_temp(FCD_IIb));
            pooled_std = sqrt(power(std_z_lesion_temp_a, 2)/num_of_FCDIIA+power(std_z_lesion_temp_b, 2)/num_of_FCDIIB);
            t = (mean_z_lesion_temp_a - mean_z_lesion_temp_b)/pooled_std;
            df = power(power(std_z_lesion_temp_a, 2)/num_of_FCDIIA+power(std_z_lesion_temp_b, 2)/num_of_FCDIIB, 2)/...
                (power(power(std_z_lesion_temp_a, 2)/num_of_FCDIIA,2)/(num_of_FCDIIA-1) + power(power(std_z_lesion_temp_b, 2)/num_of_FCDIIB,2)/(num_of_FCDIIB-1));
            p = tcdf(t, df);
            p_type_IIa_vs_IIb = [ p_type_IIa_vs_IIb min(p, 1-p) ]
%             g = [-1*ones(length(feature_temp(FCD_IIa)), 1); ones(length(feature_temp(FCD_IIb)), 1); ]; M = 1 + term(g);
%             slm = SurfStatLinMod([feature_temp(FCD_IIa); feature_temp(FCD_IIb)], M);
%             slm = SurfStatT(slm, g); p = 1 - tcdf(slm.t,slm.df);
%             p_type_IIa_vs_IIb_ss = [ p_type_IIa_vs_IIb_ss 2*min(p, 1-p) ];
            
        end
        for feature_BC = 1
            
            feature_temp = Lesion_BC;
            
            mean_z_lesion_temp = mean(feature_temp);
            std_z_lesion_temp  = std(feature_temp);
            pooled_std = sqrt(power(std_z_lesion_temp, 2)/num_of_FCD+1/num_of_cont);
            t = (mean_z_lesion_temp - 0)/pooled_std;
            df = power(power(std_z_lesion_temp, 2)/num_of_FCD+1/num_of_cont, 2)/(power(power(std_z_lesion_temp, 2)/num_of_FCD,2)/(num_of_FCD-1) + power(1/num_of_cont, 2)/(num_of_cont-1));
            p = tcdf(t, df);
            p_type = [ p_type 2*min(p, 1-p) ]
%             g = ones(length(feature_temp), 1); M = 1 + term(g);
%             slm = SurfStatLinMod(feature_temp, M);
%             slm = SurfStatT(slm, g); p = 1 - tcdf(slm.t,slm.df);
%             p_type_ss = [ p_type_ss 2*min(p, 1-p) ];
            
            mean_z_lesion_temp = mean(feature_temp(FCD_IIa));
            std_z_lesion_temp  = std(feature_temp(FCD_IIa));
            pooled_std = sqrt(power(std_z_lesion_temp, 2)/num_of_FCDIIA+1/num_of_cont);
            t = (mean_z_lesion_temp - 0)/pooled_std;
            df = power(power(std_z_lesion_temp, 2)/num_of_FCDIIA+1/num_of_cont, 2)/(power(power(std_z_lesion_temp, 2)/num_of_FCDIIA,2)/(num_of_FCDIIA-1) + power(1/num_of_cont, 2)/(num_of_cont-1));
            p = tcdf(t, df);
            p_type_IIa = [ p_type_IIa 2*min(p, 1-p) ]
%             g = ones(length(feature_temp(FCD_IIa)), 1); M = 1 + term(g);
%             slm = SurfStatLinMod(feature_temp(FCD_IIa), M);
%             slm = SurfStatT(slm, g); p = 1 - tcdf(slm.t,slm.df);
            p_type_IIa_ss = [ p_type_IIa_ss 2*min(p, 1-p) ];
            
            mean_z_lesion_temp = mean(feature_temp(FCD_IIb));
            std_z_lesion_temp  = std(feature_temp(FCD_IIb));
            pooled_std = sqrt(power(std_z_lesion_temp, 2)/num_of_FCDIIB+1/num_of_cont);
            t = (mean_z_lesion_temp - 0)/pooled_std;
            df = power(power(std_z_lesion_temp, 2)/num_of_FCDIIB+1/num_of_cont, 2)/(power(power(std_z_lesion_temp, 2)/num_of_FCDIIB,2)/(num_of_FCDIIB-1) + power(1/num_of_cont, 2)/(num_of_cont-1));
            p = tcdf(t, df);
            p_type_IIb = [ p_type_IIb 2*min(p, 1-p) ]
%             g = ones(length(feature_temp(FCD_IIb)), 1); M = 1 + term(g);
%             slm = SurfStatLinMod(feature_temp(FCD_IIb), M);
%             slm = SurfStatT(slm, g); p = 1 - tcdf(slm.t,slm.df);
%             p_type_IIb_ss = [ p_type_IIb_ss 2*min(p, 1-p) ];
            
            mean_z_lesion_temp_a = mean(feature_temp(FCD_IIa));
            std_z_lesion_temp_a  = std(feature_temp(FCD_IIa));
            mean_z_lesion_temp_b = mean(feature_temp(FCD_IIb));
            std_z_lesion_temp_b  = std(feature_temp(FCD_IIb));
            pooled_std = sqrt(power(std_z_lesion_temp_a, 2)/num_of_FCDIIA+power(std_z_lesion_temp_b, 2)/num_of_FCDIIB);
            t = (mean_z_lesion_temp_a - mean_z_lesion_temp_b)/pooled_std;
            df = power(power(std_z_lesion_temp_a, 2)/num_of_FCDIIA+power(std_z_lesion_temp_b, 2)/num_of_FCDIIB, 2)/...
                (power(power(std_z_lesion_temp_a, 2)/num_of_FCDIIA,2)/(num_of_FCDIIA-1) + power(power(std_z_lesion_temp_b, 2)/num_of_FCDIIB,2)/(num_of_FCDIIB-1));
            p = tcdf(t, df);
            p_type_IIa_vs_IIb = [ p_type_IIa_vs_IIb min(p, 1-p) ]
%             g = [-1*ones(length(feature_temp(FCD_IIa)), 1); ones(length(feature_temp(FCD_IIb)), 1); ]; M = 1 + term(g);
%             slm = SurfStatLinMod([feature_temp(FCD_IIa); feature_temp(FCD_IIb)], M);
%             slm = SurfStatT(slm, g); p = 1 - tcdf(slm.t,slm.df);
%             p_type_IIa_vs_IIb_ss = [ p_type_IIa_vs_IIb_ss 2*min(p, 1-p) ];
            
        end
        for feature_etc = 1
%             for feature_fALFF5 = 1
%                 feature_temp = Lesion_fALFF5;
%                 
%                 mean_z_lesion_temp = mean(feature_temp);
%                 std_z_lesion_temp  = std(feature_temp);
%                 pooled_std = sqrt(power(std_z_lesion_temp, 2)/num_of_FCD+1/num_of_cont);
%                 t = (mean_z_lesion_temp - 0)/pooled_std;
%                 df = power(power(std_z_lesion_temp, 2)/num_of_FCD+1/num_of_cont, 2)/(power(power(std_z_lesion_temp, 2)/num_of_FCD,2)/(num_of_FCD-1) + power(1/num_of_cont, 2)/(num_of_cont-1));
%                 p = tcdf(t, df);
%                 p_type = [ p_type 2*min(p, 1-p) ]
%                 g = ones(length(feature_temp), 1); M = 1 + term(g);
%                 slm = SurfStatLinMod(feature_temp, M);
%                 slm = SurfStatT(slm, g); p = 1 - tcdf(slm.t,slm.df);
%                 p_type_ss = [ p_type_ss 2*min(p, 1-p) ];
%                 
%                 mean_z_lesion_temp = mean(feature_temp(FCD_IIa));
%                 std_z_lesion_temp  = std(feature_temp(FCD_IIa));
%                 pooled_std = sqrt(power(std_z_lesion_temp, 2)/num_of_FCDIIA+1/num_of_cont);
%                 t = (mean_z_lesion_temp - 0)/pooled_std;
%                 df = power(power(std_z_lesion_temp, 2)/num_of_FCDIIA+1/num_of_cont, 2)/(power(power(std_z_lesion_temp, 2)/num_of_FCDIIA,2)/(num_of_FCDIIA-1) + power(1/num_of_cont, 2)/(num_of_cont-1));
%                 p = tcdf(t, df);
%                 p_type_IIa = [ p_type_IIa 2*min(p, 1-p) ]
%                 g = ones(length(feature_temp(FCD_IIa)), 1); M = 1 + term(g);
%                 slm = SurfStatLinMod(feature_temp(FCD_IIa), M);
%                 slm = SurfStatT(slm, g); p = 1 - tcdf(slm.t,slm.df);
%                 p_type_IIa_ss = [ p_type_IIa_ss 2*min(p, 1-p) ];
%                 
%                 mean_z_lesion_temp = mean(feature_temp(FCD_IIb));
%                 std_z_lesion_temp  = std(feature_temp(FCD_IIb));
%                 pooled_std = sqrt(power(std_z_lesion_temp, 2)/num_of_FCDIIB+1/num_of_cont);
%                 t = (mean_z_lesion_temp - 0)/pooled_std;
%                 df = power(power(std_z_lesion_temp, 2)/num_of_FCDIIB+1/num_of_cont, 2)/(power(power(std_z_lesion_temp, 2)/num_of_FCDIIB,2)/(num_of_FCDIIB-1) + power(1/num_of_cont, 2)/(num_of_cont-1));
%                 p = tcdf(t, df);
%                 p_type_IIb = [ p_type_IIb 2*min(p, 1-p) ]
%                 g = ones(length(feature_temp(FCD_IIb)), 1); M = 1 + term(g);
%                 slm = SurfStatLinMod(feature_temp(FCD_IIb), M);
%                 slm = SurfStatT(slm, g); p = 1 - tcdf(slm.t,slm.df);
%                 p_type_IIb_ss = [ p_type_IIb_ss 2*min(p, 1-p) ];
%                 
%                 mean_z_lesion_temp_a = mean(feature_temp(FCD_IIa));
%                 std_z_lesion_temp_a  = std(feature_temp(FCD_IIa));
%                 mean_z_lesion_temp_b = mean(feature_temp(FCD_IIb));
%                 std_z_lesion_temp_b  = std(feature_temp(FCD_IIb));
%                 pooled_std = sqrt(power(std_z_lesion_temp_a, 2)/num_of_FCDIIA+power(std_z_lesion_temp_b, 2)/num_of_FCDIIB);
%                 t = (mean_z_lesion_temp_a - mean_z_lesion_temp_b)/pooled_std;
%                 df = power(power(std_z_lesion_temp_a, 2)/num_of_FCDIIA+power(std_z_lesion_temp_b, 2)/num_of_FCDIIB, 2)/...
%                     (power(power(std_z_lesion_temp_a, 2)/num_of_FCDIIA,2)/(num_of_FCDIIA-1) + power(power(std_z_lesion_temp_b, 2)/num_of_FCDIIB,2)/(num_of_FCDIIB-1));
%                 p = tcdf(t, df);
%                 p_type_IIa_vs_IIb = [ p_type_IIa_vs_IIb min(p, 1-p) ]
%                 g = [-1*ones(length(feature_temp(FCD_IIa)), 1); ones(length(feature_temp(FCD_IIb)), 1); ]; M = 1 + term(g);
%                 slm = SurfStatLinMod([feature_temp(FCD_IIa); feature_temp(FCD_IIb)], M);
%                 slm = SurfStatT(slm, g); p = 1 - tcdf(slm.t,slm.df);
%                 p_type_IIa_vs_IIb_ss = [ p_type_IIa_vs_IIb_ss 2*min(p, 1-p) ];
%             end
%             for feature_DC5 = 1
%                 
%                 feature_temp = Lesion_DC5;
%                 
%                 mean_z_lesion_temp = mean(feature_temp);
%                 std_z_lesion_temp  = std(feature_temp);
%                 pooled_std = sqrt(power(std_z_lesion_temp, 2)/num_of_FCD+1/num_of_cont);
%                 t = (mean_z_lesion_temp - 0)/pooled_std;
%                 df = power(power(std_z_lesion_temp, 2)/num_of_FCD+1/num_of_cont, 2)/(power(power(std_z_lesion_temp, 2)/num_of_FCD,2)/(num_of_FCD-1) + power(1/num_of_cont, 2)/(num_of_cont-1));
%                 p = tcdf(t, df);
%                 p_type = [ p_type 2*min(p, 1-p) ]
%                 g = ones(length(feature_temp), 1); M = 1 + term(g);
%                 slm = SurfStatLinMod(feature_temp, M);
%                 slm = SurfStatT(slm, g); p = 1 - tcdf(slm.t,slm.df);
%                 p_type_ss = [ p_type_ss 2*min(p, 1-p) ];
%                 
%                 mean_z_lesion_temp = mean(feature_temp(FCD_IIa));
%                 std_z_lesion_temp  = std(feature_temp(FCD_IIa));
%                 pooled_std = sqrt(power(std_z_lesion_temp, 2)/num_of_FCDIIA+1/num_of_cont);
%                 t = (mean_z_lesion_temp - 0)/pooled_std;
%                 df = power(power(std_z_lesion_temp, 2)/num_of_FCDIIA+1/num_of_cont, 2)/(power(power(std_z_lesion_temp, 2)/num_of_FCDIIA,2)/(num_of_FCDIIA-1) + power(1/num_of_cont, 2)/(num_of_cont-1));
%                 p = tcdf(t, df);
%                 p_type_IIa = [ p_type_IIa 2*min(p, 1-p) ]
%                 g = ones(length(feature_temp(FCD_IIa)), 1); M = 1 + term(g);
%                 slm = SurfStatLinMod(feature_temp(FCD_IIa), M);
%                 slm = SurfStatT(slm, g); p = 1 - tcdf(slm.t,slm.df);
%                 p_type_IIa_ss = [ p_type_IIa_ss 2*min(p, 1-p) ];
%                 
%                 mean_z_lesion_temp = mean(feature_temp(FCD_IIb));
%                 std_z_lesion_temp  = std(feature_temp(FCD_IIb));
%                 pooled_std = sqrt(power(std_z_lesion_temp, 2)/num_of_FCDIIB+1/num_of_cont);
%                 t = (mean_z_lesion_temp - 0)/pooled_std;
%                 df = power(power(std_z_lesion_temp, 2)/num_of_FCDIIB+1/num_of_cont, 2)/(power(power(std_z_lesion_temp, 2)/num_of_FCDIIB,2)/(num_of_FCDIIB-1) + power(1/num_of_cont, 2)/(num_of_cont-1));
%                 p = tcdf(t, df);
%                 p_type_IIb = [ p_type_IIb 2*min(p, 1-p) ]
%                 g = ones(length(feature_temp(FCD_IIb)), 1); M = 1 + term(g);
%                 slm = SurfStatLinMod(feature_temp(FCD_IIb), M);
%                 slm = SurfStatT(slm, g); p = 1 - tcdf(slm.t,slm.df);
%                 p_type_IIb_ss = [ p_type_IIb_ss 2*min(p, 1-p) ];
%                 
%                 mean_z_lesion_temp_a = mean(feature_temp(FCD_IIa));
%                 std_z_lesion_temp_a  = std(feature_temp(FCD_IIa));
%                 mean_z_lesion_temp_b = mean(feature_temp(FCD_IIb));
%                 std_z_lesion_temp_b  = std(feature_temp(FCD_IIb));
%                 pooled_std = sqrt(power(std_z_lesion_temp_a, 2)/num_of_FCDIIA+power(std_z_lesion_temp_b, 2)/num_of_FCDIIB);
%                 t = (mean_z_lesion_temp_a - mean_z_lesion_temp_b)/pooled_std;
%                 df = power(power(std_z_lesion_temp_a, 2)/num_of_FCDIIA+power(std_z_lesion_temp_b, 2)/num_of_FCDIIB, 2)/...
%                     (power(power(std_z_lesion_temp_a, 2)/num_of_FCDIIA,2)/(num_of_FCDIIA-1) + power(power(std_z_lesion_temp_b, 2)/num_of_FCDIIB,2)/(num_of_FCDIIB-1));
%                 p = tcdf(t, df);
%                 p_type_IIa_vs_IIb = [ p_type_IIa_vs_IIb min(p, 1-p) ]
%                 g = [-1*ones(length(feature_temp(FCD_IIa)), 1); ones(length(feature_temp(FCD_IIb)), 1); ]; M = 1 + term(g);
%                 slm = SurfStatLinMod([feature_temp(FCD_IIa); feature_temp(FCD_IIb)], M);
%                 slm = SurfStatT(slm, g); p = 1 - tcdf(slm.t,slm.df);
%                 p_type_IIa_vs_IIb_ss = [ p_type_IIa_vs_IIb_ss 2*min(p, 1-p) ];
%                 
%             end
%             for feature_DC2_5 = 1
%                 
%                 feature_temp = Lesion_DC2_5;
%                 
%                 mean_z_lesion_temp = mean(feature_temp);
%                 std_z_lesion_temp  = std(feature_temp);
%                 pooled_std = sqrt(power(std_z_lesion_temp, 2)/num_of_FCD+1/num_of_cont);
%                 t = (mean_z_lesion_temp - 0)/pooled_std;
%                 df = power(power(std_z_lesion_temp, 2)/num_of_FCD+1/num_of_cont, 2)/(power(power(std_z_lesion_temp, 2)/num_of_FCD,2)/(num_of_FCD-1) + power(1/num_of_cont, 2)/(num_of_cont-1));
%                 p = tcdf(t, df);
%                 p_type = [ p_type 2*min(p, 1-p) ]
%                 g = ones(length(feature_temp), 1); M = 1 + term(g);
%                 slm = SurfStatLinMod(feature_temp, M);
%                 slm = SurfStatT(slm, g); p = 1 - tcdf(slm.t,slm.df);
%                 p_type_ss = [ p_type_ss 2*min(p, 1-p) ];
%                 
%                 mean_z_lesion_temp = mean(feature_temp(FCD_IIa));
%                 std_z_lesion_temp  = std(feature_temp(FCD_IIa));
%                 pooled_std = sqrt(power(std_z_lesion_temp, 2)/num_of_FCDIIA+1/num_of_cont);
%                 t = (mean_z_lesion_temp - 0)/pooled_std;
%                 df = power(power(std_z_lesion_temp, 2)/num_of_FCDIIA+1/num_of_cont, 2)/(power(power(std_z_lesion_temp, 2)/num_of_FCDIIA,2)/(num_of_FCDIIA-1) + power(1/num_of_cont, 2)/(num_of_cont-1));
%                 p = tcdf(t, df);
%                 p_type_IIa = [ p_type_IIa 2*min(p, 1-p) ]
%                 g = ones(length(feature_temp(FCD_IIa)), 1); M = 1 + term(g);
%                 slm = SurfStatLinMod(feature_temp(FCD_IIa), M);
%                 slm = SurfStatT(slm, g); p = 1 - tcdf(slm.t,slm.df);
%                 p_type_IIa_ss = [ p_type_IIa_ss 2*min(p, 1-p) ];
%                 
%                 mean_z_lesion_temp = mean(feature_temp(FCD_IIb));
%                 std_z_lesion_temp  = std(feature_temp(FCD_IIb));
%                 pooled_std = sqrt(power(std_z_lesion_temp, 2)/num_of_FCDIIB+1/num_of_cont);
%                 t = (mean_z_lesion_temp - 0)/pooled_std;
%                 df = power(power(std_z_lesion_temp, 2)/num_of_FCDIIB+1/num_of_cont, 2)/(power(power(std_z_lesion_temp, 2)/num_of_FCDIIB,2)/(num_of_FCDIIB-1) + power(1/num_of_cont, 2)/(num_of_cont-1));
%                 p = tcdf(t, df);
%                 p_type_IIb = [ p_type_IIb 2*min(p, 1-p) ]
%                 g = ones(length(feature_temp(FCD_IIb)), 1); M = 1 + term(g);
%                 slm = SurfStatLinMod(feature_temp(FCD_IIb), M);
%                 slm = SurfStatT(slm, g); p = 1 - tcdf(slm.t,slm.df);
%                 p_type_IIb_ss = [ p_type_IIb_ss 2*min(p, 1-p) ];
%                 
%                 mean_z_lesion_temp_a = mean(feature_temp(FCD_IIa));
%                 std_z_lesion_temp_a  = std(feature_temp(FCD_IIa));
%                 mean_z_lesion_temp_b = mean(feature_temp(FCD_IIb));
%                 std_z_lesion_temp_b  = std(feature_temp(FCD_IIb));
%                 pooled_std = sqrt(power(std_z_lesion_temp_a, 2)/num_of_FCDIIA+power(std_z_lesion_temp_b, 2)/num_of_FCDIIB);
%                 t = (mean_z_lesion_temp_a - mean_z_lesion_temp_b)/pooled_std;
%                 df = power(power(std_z_lesion_temp_a, 2)/num_of_FCDIIA+power(std_z_lesion_temp_b, 2)/num_of_FCDIIB, 2)/...
%                     (power(power(std_z_lesion_temp_a, 2)/num_of_FCDIIA,2)/(num_of_FCDIIA-1) + power(power(std_z_lesion_temp_b, 2)/num_of_FCDIIB,2)/(num_of_FCDIIB-1));
%                 p = tcdf(t, df);
%                 p_type_IIa_vs_IIb = [ p_type_IIa_vs_IIb min(p, 1-p) ]
%                 g = [-1*ones(length(feature_temp(FCD_IIa)), 1); ones(length(feature_temp(FCD_IIb)), 1); ]; M = 1 + term(g);
%                 slm = SurfStatLinMod([feature_temp(FCD_IIa); feature_temp(FCD_IIb)], M);
%                 slm = SurfStatT(slm, g); p = 1 - tcdf(slm.t,slm.df);
%                 p_type_IIa_vs_IIb_ss = [ p_type_IIa_vs_IIb_ss 2*min(p, 1-p) ];
%                 
%             end
%             for feature_DC3_5 = 1
%                 
%                 feature_temp = Lesion_DC3_5;
%                 
%                 mean_z_lesion_temp = mean(feature_temp);
%                 std_z_lesion_temp  = std(feature_temp);
%                 pooled_std = sqrt(power(std_z_lesion_temp, 2)/num_of_FCD+1/num_of_cont);
%                 t = (mean_z_lesion_temp - 0)/pooled_std;
%                 df = power(power(std_z_lesion_temp, 2)/num_of_FCD+1/num_of_cont, 2)/(power(power(std_z_lesion_temp, 2)/num_of_FCD,2)/(num_of_FCD-1) + power(1/num_of_cont, 2)/(num_of_cont-1));
%                 p = tcdf(t, df);
%                 p_type = [ p_type 2*min(p, 1-p) ]
%                 g = ones(length(feature_temp), 1); M = 1 + term(g);
%                 slm = SurfStatLinMod(feature_temp, M);
%                 slm = SurfStatT(slm, g); p = 1 - tcdf(slm.t,slm.df);
%                 p_type_ss = [ p_type_ss 2*min(p, 1-p) ];
%                 
%                 mean_z_lesion_temp = mean(feature_temp(FCD_IIa));
%                 std_z_lesion_temp  = std(feature_temp(FCD_IIa));
%                 pooled_std = sqrt(power(std_z_lesion_temp, 2)/num_of_FCDIIA+1/num_of_cont);
%                 t = (mean_z_lesion_temp - 0)/pooled_std;
%                 df = power(power(std_z_lesion_temp, 2)/num_of_FCDIIA+1/num_of_cont, 2)/(power(power(std_z_lesion_temp, 2)/num_of_FCDIIA,2)/(num_of_FCDIIA-1) + power(1/num_of_cont, 2)/(num_of_cont-1));
%                 p = tcdf(t, df);
%                 p_type_IIa = [ p_type_IIa 2*min(p, 1-p) ]
%                 g = ones(length(feature_temp(FCD_IIa)), 1); M = 1 + term(g);
%                 slm = SurfStatLinMod(feature_temp(FCD_IIa), M);
%                 slm = SurfStatT(slm, g); p = 1 - tcdf(slm.t,slm.df);
%                 p_type_IIa_ss = [ p_type_IIa_ss 2*min(p, 1-p) ];
%                 
%                 mean_z_lesion_temp = mean(feature_temp(FCD_IIb));
%                 std_z_lesion_temp  = std(feature_temp(FCD_IIb));
%                 pooled_std = sqrt(power(std_z_lesion_temp, 2)/num_of_FCDIIB+1/num_of_cont);
%                 t = (mean_z_lesion_temp - 0)/pooled_std;
%                 df = power(power(std_z_lesion_temp, 2)/num_of_FCDIIB+1/num_of_cont, 2)/(power(power(std_z_lesion_temp, 2)/num_of_FCDIIB,2)/(num_of_FCDIIB-1) + power(1/num_of_cont, 2)/(num_of_cont-1));
%                 p = tcdf(t, df);
%                 p_type_IIb = [ p_type_IIb 2*min(p, 1-p) ]
%                 g = ones(length(feature_temp(FCD_IIb)), 1); M = 1 + term(g);
%                 slm = SurfStatLinMod(feature_temp(FCD_IIb), M);
%                 slm = SurfStatT(slm, g); p = 1 - tcdf(slm.t,slm.df);
%                 p_type_IIb_ss = [ p_type_IIb_ss 2*min(p, 1-p) ];
%                 
%                 mean_z_lesion_temp_a = mean(feature_temp(FCD_IIa));
%                 std_z_lesion_temp_a  = std(feature_temp(FCD_IIa));
%                 mean_z_lesion_temp_b = mean(feature_temp(FCD_IIb));
%                 std_z_lesion_temp_b  = std(feature_temp(FCD_IIb));
%                 pooled_std = sqrt(power(std_z_lesion_temp_a, 2)/num_of_FCDIIA+power(std_z_lesion_temp_b, 2)/num_of_FCDIIB);
%                 t = (mean_z_lesion_temp_a - mean_z_lesion_temp_b)/pooled_std;
%                 df = power(power(std_z_lesion_temp_a, 2)/num_of_FCDIIA+power(std_z_lesion_temp_b, 2)/num_of_FCDIIB, 2)/...
%                     (power(power(std_z_lesion_temp_a, 2)/num_of_FCDIIA,2)/(num_of_FCDIIA-1) + power(power(std_z_lesion_temp_b, 2)/num_of_FCDIIB,2)/(num_of_FCDIIB-1));
%                 p = tcdf(t, df);
%                 p_type_IIa_vs_IIb = [ p_type_IIa_vs_IIb min(p, 1-p) ]
%                 g = [-1*ones(length(feature_temp(FCD_IIa)), 1); ones(length(feature_temp(FCD_IIb)), 1); ]; M = 1 + term(g);
%                 slm = SurfStatLinMod([feature_temp(FCD_IIa); feature_temp(FCD_IIb)], M);
%                 slm = SurfStatT(slm, g); p = 1 - tcdf(slm.t,slm.df);
%                 p_type_IIa_vs_IIb_ss = [ p_type_IIa_vs_IIb_ss 2*min(p, 1-p) ];
%                 
%             end
%             for feature_EC5 = 1
%                 
%                 feature_temp = Lesion_EC5;
%                 
%                 mean_z_lesion_temp = mean(feature_temp);
%                 std_z_lesion_temp  = std(feature_temp);
%                 pooled_std = sqrt(power(std_z_lesion_temp, 2)/num_of_FCD+1/num_of_cont);
%                 t = (mean_z_lesion_temp - 0)/pooled_std;
%                 df = power(power(std_z_lesion_temp, 2)/num_of_FCD+1/num_of_cont, 2)/(power(power(std_z_lesion_temp, 2)/num_of_FCD,2)/(num_of_FCD-1) + power(1/num_of_cont, 2)/(num_of_cont-1));
%                 p = tcdf(t, df);
%                 p_type = [ p_type 2*min(p, 1-p) ]
%                 g = ones(length(feature_temp), 1); M = 1 + term(g);
%                 slm = SurfStatLinMod(feature_temp, M);
%                 slm = SurfStatT(slm, g); p = 1 - tcdf(slm.t,slm.df);
%                 p_type_ss = [ p_type_ss 2*min(p, 1-p) ];
%                 
%                 mean_z_lesion_temp = mean(feature_temp(FCD_IIa));
%                 std_z_lesion_temp  = std(feature_temp(FCD_IIa));
%                 pooled_std = sqrt(power(std_z_lesion_temp, 2)/num_of_FCDIIA+1/num_of_cont);
%                 t = (mean_z_lesion_temp - 0)/pooled_std;
%                 df = power(power(std_z_lesion_temp, 2)/num_of_FCDIIA+1/num_of_cont, 2)/(power(power(std_z_lesion_temp, 2)/num_of_FCDIIA,2)/(num_of_FCDIIA-1) + power(1/num_of_cont, 2)/(num_of_cont-1));
%                 p = tcdf(t, df);
%                 p_type_IIa = [ p_type_IIa 2*min(p, 1-p) ]
%                 g = ones(length(feature_temp(FCD_IIa)), 1); M = 1 + term(g);
%                 slm = SurfStatLinMod(feature_temp(FCD_IIa), M);
%                 slm = SurfStatT(slm, g); p = 1 - tcdf(slm.t,slm.df);
%                 p_type_IIa_ss = [ p_type_IIa_ss 2*min(p, 1-p) ];
%                 
%                 mean_z_lesion_temp = mean(feature_temp(FCD_IIb));
%                 std_z_lesion_temp  = std(feature_temp(FCD_IIb));
%                 pooled_std = sqrt(power(std_z_lesion_temp, 2)/num_of_FCDIIB+1/num_of_cont);
%                 t = (mean_z_lesion_temp - 0)/pooled_std;
%                 df = power(power(std_z_lesion_temp, 2)/num_of_FCDIIB+1/num_of_cont, 2)/(power(power(std_z_lesion_temp, 2)/num_of_FCDIIB,2)/(num_of_FCDIIB-1) + power(1/num_of_cont, 2)/(num_of_cont-1));
%                 p = tcdf(t, df);
%                 p_type_IIb = [ p_type_IIb 2*min(p, 1-p) ]
%                 g = ones(length(feature_temp(FCD_IIb)), 1); M = 1 + term(g);
%                 slm = SurfStatLinMod(feature_temp(FCD_IIb), M);
%                 slm = SurfStatT(slm, g); p = 1 - tcdf(slm.t,slm.df);
%                 p_type_IIb_ss = [ p_type_IIb_ss 2*min(p, 1-p) ];
%                 
%                 mean_z_lesion_temp_a = mean(feature_temp(FCD_IIa));
%                 std_z_lesion_temp_a  = std(feature_temp(FCD_IIa));
%                 mean_z_lesion_temp_b = mean(feature_temp(FCD_IIb));
%                 std_z_lesion_temp_b  = std(feature_temp(FCD_IIb));
%                 pooled_std = sqrt(power(std_z_lesion_temp_a, 2)/num_of_FCDIIA+power(std_z_lesion_temp_b, 2)/num_of_FCDIIB);
%                 t = (mean_z_lesion_temp_a - mean_z_lesion_temp_b)/pooled_std;
%                 df = power(power(std_z_lesion_temp_a, 2)/num_of_FCDIIA+power(std_z_lesion_temp_b, 2)/num_of_FCDIIB, 2)/...
%                     (power(power(std_z_lesion_temp_a, 2)/num_of_FCDIIA,2)/(num_of_FCDIIA-1) + power(power(std_z_lesion_temp_b, 2)/num_of_FCDIIB,2)/(num_of_FCDIIB-1));
%                 p = tcdf(t, df);
%                 p_type_IIa_vs_IIb = [ p_type_IIa_vs_IIb min(p, 1-p) ]
%                 g = [-1*ones(length(feature_temp(FCD_IIa)), 1); ones(length(feature_temp(FCD_IIb)), 1); ]; M = 1 + term(g);
%                 slm = SurfStatLinMod([feature_temp(FCD_IIa); feature_temp(FCD_IIb)], M);
%                 slm = SurfStatT(slm, g); p = 1 - tcdf(slm.t,slm.df);
%                 p_type_IIa_vs_IIb_ss = [ p_type_IIa_vs_IIb_ss 2*min(p, 1-p) ];
%                 
%             end
%             for feature_VMHC5 = 1
%                 
%                 feature_temp = Lesion_VMHC5;
%                 
%                 mean_z_lesion_temp = mean(feature_temp);
%                 std_z_lesion_temp  = std(feature_temp);
%                 pooled_std = sqrt(power(std_z_lesion_temp, 2)/num_of_FCD+1/num_of_cont);
%                 t = (mean_z_lesion_temp - 0)/pooled_std;
%                 df = power(power(std_z_lesion_temp, 2)/num_of_FCD+1/num_of_cont, 2)/(power(power(std_z_lesion_temp, 2)/num_of_FCD,2)/(num_of_FCD-1) + power(1/num_of_cont, 2)/(num_of_cont-1));
%                 p = tcdf(t, df);
%                 p_type = [ p_type 2*min(p, 1-p) ]
%                 g = ones(length(feature_temp), 1); M = 1 + term(g);
%                 slm = SurfStatLinMod(feature_temp, M);
%                 slm = SurfStatT(slm, g); p = 1 - tcdf(slm.t,slm.df);
%                 p_type_ss = [ p_type_ss 2*min(p, 1-p) ];
%                 
%                 mean_z_lesion_temp = mean(feature_temp(FCD_IIa));
%                 std_z_lesion_temp  = std(feature_temp(FCD_IIa));
%                 pooled_std = sqrt(power(std_z_lesion_temp, 2)/num_of_FCDIIA+1/num_of_cont);
%                 t = (mean_z_lesion_temp - 0)/pooled_std;
%                 df = power(power(std_z_lesion_temp, 2)/num_of_FCDIIA+1/num_of_cont, 2)/(power(power(std_z_lesion_temp, 2)/num_of_FCDIIA,2)/(num_of_FCDIIA-1) + power(1/num_of_cont, 2)/(num_of_cont-1));
%                 p = tcdf(t, df);
%                 p_type_IIa = [ p_type_IIa 2*min(p, 1-p) ]
%                 g = ones(length(feature_temp(FCD_IIa)), 1); M = 1 + term(g);
%                 slm = SurfStatLinMod(feature_temp(FCD_IIa), M);
%                 slm = SurfStatT(slm, g); p = 1 - tcdf(slm.t,slm.df);
%                 p_type_IIa_ss = [ p_type_IIa_ss 2*min(p, 1-p) ];
%                 
%                 mean_z_lesion_temp = mean(feature_temp(FCD_IIb));
%                 std_z_lesion_temp  = std(feature_temp(FCD_IIb));
%                 pooled_std = sqrt(power(std_z_lesion_temp, 2)/num_of_FCDIIB+1/num_of_cont);
%                 t = (mean_z_lesion_temp - 0)/pooled_std;
%                 df = power(power(std_z_lesion_temp, 2)/num_of_FCDIIB+1/num_of_cont, 2)/(power(power(std_z_lesion_temp, 2)/num_of_FCDIIB,2)/(num_of_FCDIIB-1) + power(1/num_of_cont, 2)/(num_of_cont-1));
%                 p = tcdf(t, df);
%                 p_type_IIb = [ p_type_IIb 2*min(p, 1-p) ]
%                 g = ones(length(feature_temp(FCD_IIb)), 1); M = 1 + term(g);
%                 slm = SurfStatLinMod(feature_temp(FCD_IIb), M);
%                 slm = SurfStatT(slm, g); p = 1 - tcdf(slm.t,slm.df);
%                 p_type_IIb_ss = [ p_type_IIb_ss 2*min(p, 1-p) ];
%                 
%                 mean_z_lesion_temp_a = mean(feature_temp(FCD_IIa));
%                 std_z_lesion_temp_a  = std(feature_temp(FCD_IIa));
%                 mean_z_lesion_temp_b = mean(feature_temp(FCD_IIb));
%                 std_z_lesion_temp_b  = std(feature_temp(FCD_IIb));
%                 pooled_std = sqrt(power(std_z_lesion_temp_a, 2)/num_of_FCDIIA+power(std_z_lesion_temp_b, 2)/num_of_FCDIIB);
%                 t = (mean_z_lesion_temp_a - mean_z_lesion_temp_b)/pooled_std;
%                 df = power(power(std_z_lesion_temp_a, 2)/num_of_FCDIIA+power(std_z_lesion_temp_b, 2)/num_of_FCDIIB, 2)/...
%                     (power(power(std_z_lesion_temp_a, 2)/num_of_FCDIIA,2)/(num_of_FCDIIA-1) + power(power(std_z_lesion_temp_b, 2)/num_of_FCDIIB,2)/(num_of_FCDIIB-1));
%                 p = tcdf(t, df);
%                 p_type_IIa_vs_IIb = [ p_type_IIa_vs_IIb min(p, 1-p) ]
%                 g = [-1*ones(length(feature_temp(FCD_IIa)), 1); ones(length(feature_temp(FCD_IIb)), 1); ]; M = 1 + term(g);
%                 slm = SurfStatLinMod([feature_temp(FCD_IIa); feature_temp(FCD_IIb)], M);
%                 slm = SurfStatT(slm, g); p = 1 - tcdf(slm.t,slm.df);
%                 p_type_IIa_vs_IIb_ss = [ p_type_IIa_vs_IIb_ss 2*min(p, 1-p) ];
%                 
            end
    end
    
%     FDR_threshold_final_lesion_profiling = FDR([p_value_lesion_profile_T1_FLAIR p_value_lesion_profile_DTI p_type_IIa p_type_IIb p_type_IIa_vs_IIb], 0.05); %% 0.0128
    FDR_threshold_final_lesion_profiling = FDR(p_value_set_all_lesion_profiling, 0.05); %% 0.02
    FDR_threshold = FDR_threshold_final_lesion_profiling;
    
    for visualization = 1
        
        for depth_format = 1
                        
            feature_profile_typeIIa = [                
                Lesion_ALFF5(FCD_IIa)'
                Lesion_ReHo5(FCD_IIa)'
                Lesion_DC(FCD_IIa)
                Lesion_BC(FCD_IIa)                
                ]';
            
            feature_profile_typeIIb = [                
                Lesion_ALFF5(FCD_IIb)'
                Lesion_ReHo5(FCD_IIb)'
                Lesion_DC(FCD_IIb)
                Lesion_BC(FCD_IIb)                
                ]';
            
            feature_idx = [ 1 ];
            x_lim = [0 3+(length(feature_idx)-1)*4];
            y_lim = [-2 2];
            p_value_set_ALFF = [ p_type_IIa(feature_idx) p_type_IIb(feature_idx) ];
            visualization_graph_cortical_depth_profile(feature_idx, p_value_set_ALFF, feature_profile_typeIIa(:, feature_idx), feature_profile_typeIIb(:, feature_idx), x_lim, y_lim, FDR_threshold);
            
            feature_idx = [ 2 ];
            x_lim = [0 3+(length(feature_idx)-1)*4];
            y_lim = [-3 3];
            p_value_set_ReHo = [ p_type_IIa(feature_idx) p_type_IIb(feature_idx) ];
            visualization_graph_cortical_depth_profile(feature_idx, p_value_set_ReHo, feature_profile_typeIIa(:, feature_idx), feature_profile_typeIIb(:, feature_idx), x_lim, y_lim, FDR_threshold);
            
        end
        
        for global_format = 1
            
            feature_profile_typeIIa = [
                ones(1, size(find(FCD_IIa), 2))*NaN
                Lesion_ALFF5(FCD_IIa)'
                Lesion_ReHo5(FCD_IIa)'
                Lesion_DC(FCD_IIa)
                Lesion_BC(FCD_IIa)
                ones(1, size(find(FCD_IIa), 2))*NaN
                ]';
            
            feature_profile_typeIIb = [
                ones(1, size(find(FCD_IIb), 2))*NaN
                Lesion_ALFF5(FCD_IIb)'
                Lesion_ReHo5(FCD_IIb)'
                Lesion_DC(FCD_IIb)
                Lesion_BC(FCD_IIb)
                ones(1, size(find(FCD_IIb), 2))*NaN
                ]';
            
            f = figure; set(gcf, 'Position', [587 500 600 500]);
            hold on;
            linecolor = [0.8 0.8 0.8];
            hp1 = plot([0 16], [0 0]);
            hp2 = plot([0 16], [2 2], '--');
            hp3 = plot([0 16], [-2 -2], '--');
            set(hp1, 'Color', linecolor);
            set(hp2, 'Color', linecolor);
            set(hp3, 'Color', linecolor);
            set(gcf, 'Color', 'w');
            xtick = [0 3 6 9 12 16];
            h1 = notBoxPlot(feature_profile_typeIIa, xtick, 0.5, 'line');
            h2 = notBoxPlot(feature_profile_typeIIb, xtick+1, 0.5, 'line');
            flag_mu  = 'off';
            flag_sd  = 'off';
            flag_sem = 'off';
            barhead_linewidth = 0.3;
            for i = 2 : size(h1, 2)
                set(h1(i).data, 'markerfacecolor',[1,1,1],'color',[0,0,0]);
                set(h1(i).mu,   'visible' , flag_mu, 'markerfacecolor',[0,0,0],'color',[0,0,0]);
                set(h1(i).sd,   'visible' , flag_sd, 'LineStyle', ':', 'Color',[0.3,0.3,0.3]);
                set(h1(i).sem,  'visible' , flag_sem, 'LineStyle', '-', 'Color',[0.3,0.3,0.3]);
                set(h2(i).mu,   'visible' , flag_mu, 'markerfacecolor',[0,0,0],'color',[0,0,0]);
                set(h2(i).sd,   'visible' , flag_sd, 'LineStyle', ':', 'Color',[0.3,0.3,0.3]);
                set(h2(i).sem,  'visible' , flag_sem, 'LineStyle', '-', 'Color',[0.3,0.3,0.3]);
                
                x_1_temp = get(h1(i).mu, 'xData');
                x_2_temp = get(h2(i).mu, 'xData');
                mu_1 = get(h1(i).mu, 'YData');
                mu_2 = get(h2(i).mu, 'YData');
                sd_1 = get(h1(i).sd, 'YData');
                sd_2 = get(h2(i).sd, 'YData');
                med_1 = median(get(h1(i).data, 'YData'));
                med_2 = median(get(h2(i).data, 'YData'));
                
                hold on; line([x_1_temp-barhead_linewidth x_1_temp+barhead_linewidth], [sd_1(1) sd_1(1)], 'LineWidth',2,'Color',[0 0 0]);
                hold on; line([x_1_temp-barhead_linewidth x_1_temp+barhead_linewidth], [sd_1(2) sd_1(2)], 'LineWidth',2,'Color',[0 0 0]);
                hold on; line([x_2_temp-barhead_linewidth x_2_temp+barhead_linewidth], [sd_2(1) sd_2(1)], 'LineWidth',2,'Color',[0 0 0]);
                hold on; line([x_2_temp-barhead_linewidth x_2_temp+barhead_linewidth], [sd_2(2) sd_2(2)], 'LineWidth',2,'Color',[0 0 0]);
                hold on; line([x_1_temp-barhead_linewidth x_1_temp+barhead_linewidth], [mu_1 mu_1], 'LineWidth',2,'Color',[0 0 0]);
                hold on; line([x_2_temp-barhead_linewidth x_2_temp+barhead_linewidth], [mu_2 mu_2], 'LineWidth',2,'Color',[0 0 0]);
                hold on; scatter(x_1_temp, mu_1, 60, 'MarkerEdgeColor', [0 0 0], 'markerfacecolor', [0 0 0]);
                hold on; scatter(x_2_temp, mu_2, 60, 'MarkerEdgeColor', [0 0 0], 'markerfacecolor', [0 0 0]);
            end
            
            xlim([2 max(xtick)-2]);
            set(gca, 'XTick', xtick+0.5);
            label_feature_profile_type = {'', 'ALFF', 'ReHo', 'DegCen', 'BetCen', ''};
            set(gca, 'XTickLabel', label_feature_profile_type);
            ylabel('z scores');
            ylim([-5 5]);
            set(gca, 'YTick', [-5 -2 0 2 5]);
            
            x_range = xtick(2:end-1);
            FDR_threshold = FDR_threshold_final_lesion_profiling;
            for i = 1 : size(p_type_IIa, 2)
                if(p_type_IIa(i) <= FDR_threshold)
                    scatter(x_range(i), 3, 50, 'rv');
                elseif(p_type_IIa(i) < 0.05 && p_type_IIa(i) > FDR_threshold)
                    scatter(x_range(i), 3, 50, 'r*');
                end
            end
            
            for i = 1 : size(p_type_IIb, 2)
                if(p_type_IIb(i) <= FDR_threshold)
                    scatter(x_range(i)+1, 3, 50, 'rv');
                elseif(p_type_IIb(i) < 0.05 && p_type_IIb(i) > FDR_threshold)
                    scatter(x_range(i)+1, 3, 50, 'r*');
                end
            end
            
        end
        
        if(save_file)
            export_fig([OUTPATH '/38_FCD_lesion_feature_profile_typeIIa_IIb_ALFF' ], '-m4', '-png'); close(gcf);
        end
        
    end
    
end