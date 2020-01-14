clear all
close all
addpath('/data/noel/noel2/local/matlab/surfstat_chicago');

%% The difference with the original script is that here, we have control data that were averaged across lesion locations.
%% Therefore we are able to run non-parametric tests. But you should note that this way to generate control data can artificially
%% amplify the statistical significances as it actually reduces the standard deviation of controls to less than 1.

%% Read pre-calculated zscore of T1 and FLAIR in the lesion
for read_matfiles = 1
    
    load('zscore_database_T1_FLAIR.mat');
    load('zscore_database_T1_FLAIR_cont.mat');
    load('p_value_set_lesion_profile_T1_FLAIR.mat');
    load('p_value_set_lesion_profile_DTI.mat');
    load('p_value_set_lesion_cortical_depth_profile_T1_FLAIR.mat');
    load('p_value_set_lesion_cortical_depth_profile_DTI.mat');
    load('/local_raid/seokjun/01_project/04_IntracorticalAnalysis/01_analysis/pipeline/postprocessing_zscore_visualization/colormap_noel1.mat');
    
end

%% Read demograpic data and set up some parameters and variables
for read_demodata_setup_params = 1
    
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
    Kernel                    = 5
    Parametric                = 'quadratic'
    average_surface_dir       = '/local_raid/seokjun/01_project/04_IntracorticalAnalysis/03_result/average_surfaces/'
    
    visualization = 1;
    save_file = 0;
    visualization_linearmodel = 0;
    save_file_linearmodel = 0;
    
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
    
    % Age test: FCD vs. Controls
    [h,p,ci,stats] = ttest2(age_cont,age_pat)
    [h,p,ci,stats] = ttest2(age_cont,age_pat(strcmp(histo_type, 'FCDIIa')))
    [h,p,ci,stats] = ttest2(age_cont,age_pat(strcmp(histo_type, 'FCDIIb')))
    
    % gender test: FCD vs. Controls
    sex_tab = [ sum(strcmp(gender_cont, 'f')) sum(strcmp(gender_cont, 'm'));
                sum(strcmp(gender_pat, 'f')) sum(strcmp(gender_pat, 'm')); ];
    [p,x2] = chisquarecont(sex_tab)
    
    sex_tab = [ sum(strcmp(gender_cont, 'f')) sum(strcmp(gender_cont, 'm'));
                sum(strcmp(gender_pat(strcmp(histo_type, 'FCDIIa')), 'f')) sum(strcmp(gender_pat(strcmp(histo_type, 'FCDIIa')), 'm')); ];
    [p,x2] = chisquarecont(sex_tab)
    
    sex_tab = [ sum(strcmp(gender_cont, 'f')) sum(strcmp(gender_cont, 'm'));
                sum(strcmp(gender_pat(strcmp(histo_type, 'FCDIIb')), 'f')) sum(strcmp(gender_pat(strcmp(histo_type, 'FCDIIb')), 'm')); ];
    [p,x2] = chisquarecont(sex_tab)
    
    % Age test: FCD Type IIa vs. Type IIb
    [h,p,ci,stats] = ttest2(age_pat(strcmp(histo_type, 'FCDIIa')), age_pat(strcmp(histo_type, 'FCDIIb')))
    
    % gender test: FCD Type IIa vs. Type IIb
    sex_tab = [ sum(strcmp(gender_pat(strcmp(histo_type, 'FCDIIa')), 'f')) sum(strcmp(gender_pat(strcmp(histo_type, 'FCDIIa')), 'm'));
        sum(strcmp(gender_pat(strcmp(histo_type, 'FCDIIb')), 'f')) sum(strcmp(gender_pat(strcmp(histo_type, 'FCDIIb')), 'm')); ];
    [p,x2] = chisquarecont(sex_tab)
    
        
    FDR_threshold_final_lesion_profiling = FDR([p_value_lesion_profile_T1_FLAIR p_value_lesion_profile_DTI], 0.05); %% 0.0114
    FDR_threshold_final_cortical_depth_profiling = FDR([p_value_lesion_cortical_depth_profile p_value_lesion_cortical_depth_profile_DTI], 0.05); %% 0.02 
    
    OUTPATH = '/local_raid/seokjun/01_project/04_IntracorticalAnalysis/03_result/lesion_profile/';
    
end

%% Estimate the volume threshold to separate lesions into small, subtle and large, visible ones
for volume_thres_lesion = 1
    
    %% Using entropy index, calculate the threshold for separating lesions into small and large groups
    [a b] = sort(lesion_volume);
    Gain       = zeros(1, length(lesion_volume));
    entropy_sm = zeros(1, length(lesion_volume));
    entropy_lg = zeros(1, length(lesion_volume));
    count = 1;
    initial_temp = initial;
    
    %% One is assigned to all the volume greater than 1366
    %% because even though using original volume, 1366 is the volume giving
    %% same maximum information gain.
    initial(lesion_volume>1366) = 1;
    
    for i = 1:length(a)
        small_group = find(lesion_volume <= a(i));
        large_group = find(lesion_volume >  a(i));
        
        p = sum(initial(small_group))/length(small_group);
        if(p ~= 0 && p ~= 1 && ~isnan(p))
            entropy_sm(count) = -p*log2(p) - (1 - p)*log2(1 - p);
        else
            entropy_sm(count) = 0;
        end
        p = sum(initial(large_group))/length(large_group);
        if(p ~= 0 && p ~= 1 && ~isnan(p))
            entropy_lg(count) = -p*log2(p) - (1 - p)*log2(1 - p);
        else
            entropy_lg(count) = 0;
        end
        p = sum(initial)/length(initial);
        entropy = -p*log2(p) - (1 - p)*log2(1 - p);
        
        Gain(count) = entropy - (entropy_sm(count) + entropy_lg(count));
        count = count + 1;
    end
    initial = initial_temp;
    [C I] = max(Gain);
    
    if(visualization)
        f = figure; plot(a, Gain); xlim([0 10000]);
        hold on; plot([0 max(a)], [0 0], 'k'); plot([a(I) a(I)], [0 C], 'r--');
        s = scatter(double(a(I)), C, 'ro', 'filled'); set(s, 'SizeData', 12)
        xlabel('Threshold position along volume continuum (cm3)');
        ylabel('Information Gain');
        set(gca, 'XTickLabel', {'0', '', '2', '', '4', '', '6', '', '8', '', '10'});
        text(double(a(I))+100, C, ['T=' num2str(a(I)) 'mm^3']);
        if(save_file)
            export_fig([OUTPATH '/0_Information_Gain_Threshold' ], '-m4', '-png'); close(gcf);
        end
    end
    
    volume_threshold = a(I);
    
end

%% Select features to be included in the analyses
for feature_selection = 1
    
    %% mean_z_lesion
    %% 1:5   = RI of T1                -> 2:5      p1
    %% 6:10  = PG of T1                -> 7:9 16   p2
    %% 11:15 = TG of T1                -> 12:15    p3
    %% 16    = PG_gw of T1             -> 7:9 16   p2
    %% 17    = CT of T1                -> 17       p7
    %% 18    = MC of T1                -> 18       p8
    %% 19    = SD of T1                -> 19       p9
    %% 20:22 = RI of T1 subcortical    -> 20:22    p4
    %% 23:25 = PG of T1 subcortical    -> 23:25    p5
    %% 26:28 = TG of T1 subcortical    -> 26:28    p6
    %% 29:33 = RI of FLAIR             -> 30:33    p10
    %% 34:38 = PG of FLAIR             -> 35:37 44 p11
    %% 39:43 = TG of FLAIR             -> 40:43    p12
    %% 44    = PG_gw of FLAIR          -> 35:37 44 p11
    %% 45:47 = RI of FLAIR subcortical -> 45:47    p13
    %% 48:50 = PG of FLAIR subcortical -> 48:50    p14
    %% 51:53 = TG of FLAIR subcortical -> 51:53    p15
    feature_set_idx = { 2:5, [7:9 16], 12:15, 20:22, 17, 18, 19, 30:33, [35:37 44], 40:43, 45:47 };
    include = ones(size(case_num_pat_temp, 1), 1);
    
end

%% Start profile analyses
for profile_analyses = 1
    
    for patient_profling = 1
        
        % 1.1 Whole group patterns
        for all_patients_analysis = 1
            
            p_type_II = [];
            for i = 1 : length(feature_set_idx)
                
                pat_data = mean(mean_z_lesion(feature_set_idx{i},:), 1);
                cont_data_temp = reshape(mean(mean_z_lesion_cont(feature_set_idx{i},:, :), 2), size(mean(mean_z_lesion_cont(feature_set_idx{i},:, :), 2), 1), size(mean(mean_z_lesion_cont(feature_set_idx{i},:, :), 2), 3));
                cont_data = mean(cont_data_temp, 1);
                
                [ p h stats ] = ranksum(pat_data, cont_data);
                
                p_type_II = [ p_type_II p ];
                
            end
            p_value_set_0 = p_type_II;
            
            feature_profile = [
                ones(1, length(case_num_pat_temp))*NaN
                mean(mean_z_lesion(2:5, find(include ~= 0)), 1)
                mean(mean_z_lesion([7:9 16], find(include ~= 0)), 1)
                mean(mean_z_lesion(12:15, find(include ~= 0)), 1)
                mean(mean_z_lesion(20:22, find(include ~= 0)), 1)
                
                mean(mean_z_lesion(17, find(include ~= 0)), 1)
                mean(mean_z_lesion(18, find(include ~= 0)), 1)
                mean(mean_z_lesion(19, find(include ~= 0)), 1)
                
                mean(mean_z_lesion(30:33, find(include ~= 0)), 1)
                mean(mean_z_lesion([35:37 44], find(include ~= 0)), 1)
                mean(mean_z_lesion(40:43, find(include ~= 0)), 1)
                mean(mean_z_lesion(45:47, find(include ~= 0)), 1)
                ]';
            feature_profile(feature_profile > 100) = 10;
            
            label_feature_profile = {'Int', 'P-grad', 'T-grad', 'Sub Int', '', 'Thick', 'Curv', 'Depth', '', '', '', 'Int', 'P-grad', 'T-grad', 'Sub Int'};
            
            if(visualization)
                figure; set(gcf, 'Position', [587 500 1200 500]);
                hold on;
                linecolor = [0.8 0.8 0.8];
                hp1 = plot([0 52], [0 0]);
                hp2 = plot([0 52], [2 2], '--');
                hp3 = plot([0 52], [-2 -2], '--');
                set(hp1, 'Color', linecolor);
                set(hp2, 'Color', linecolor);
                set(hp3, 'Color', linecolor);
                set(gcf, 'Color', 'w');
                xtick = [0 3 6 9 12 16 19 22 27 30 33 36];
                h1 = notBoxPlot(feature_profile, xtick, 0.5, 'line');
                flag_mu  = 'off';
                flag_sd  = 'off';
                flag_sem = 'off';
                barhead_linewidth = 0.3;
                for i = 2 : size(h1, 2)
                    set(h1(i).data, 'markerfacecolor',[1,1,1],'color',[0,0,0]);
                    set(h1(i).mu,   'visible' , flag_mu, 'markerfacecolor',[0,0,0],'color',[0,0,0]);
                    set(h1(i).sd,   'visible' , flag_sd, 'LineStyle', ':', 'Color',[0.3,0.3,0.3]);
                    set(h1(i).sem,  'visible' , flag_sem, 'LineStyle', '-', 'Color',[0.3,0.3,0.3]);
                    
                    x_1_temp = get(h1(i).mu, 'xData');
                    mu_1 = get(h1(i).mu, 'YData');
                    sd_1 = get(h1(i).sd, 'YData');
                    med_1 = median(get(h1(i).data, 'YData'));
                    
                    hold on; line([x_1_temp-barhead_linewidth x_1_temp+barhead_linewidth], [sd_1(1) sd_1(1)], 'LineWidth',2,'Color',[0 0 0]);
                    hold on; line([x_1_temp-barhead_linewidth x_1_temp+barhead_linewidth], [sd_1(2) sd_1(2)], 'LineWidth',2,'Color',[0 0 0]);
                    hold on; line([x_1_temp-barhead_linewidth x_1_temp+barhead_linewidth], [mu_1 mu_1], 'LineWidth',2,'Color',[0 0 0]);
                    hold on; scatter(x_1_temp, mu_1, 60, 'MarkerEdgeColor', [0 0 0], 'markerfacecolor', [0 0 0]);
                end
                
                xlim([2 max(xtick)+2]);
                set(gca, 'XTick', xtick);
                %% 0    3       6         9         12        16       19      22       27      30        33        36
                label_feature_profile_type = {'', 'Int', 'P-grad', 'T-grad', 'Sub Int', 'Thick', 'Curv', 'Depth', 'Int', 'P-grad', 'T-grad', 'Sub Int'};
                set(gca, 'XTickLabel', label_feature_profile_type);
                rotateticklabel(gca, 90);
                ylabel('z scores');
                ylim([-5 10]);
                set(gca, 'YTick', [-5 -2 0 2 5]);
                
                hold on;
                x_range = xtick(2:end);
                FDR_threshold = FDR_threshold_final_lesion_profiling;
                for i = 1 : size(p_value_set_0, 2)
                    if(p_value_set_0(i) <= FDR_threshold)
                        scatter(x_range(i), 9, 50, 'rv');
                    elseif(p_value_set_0(i) < 0.05 && p_value_set_0(i) > FDR_threshold)
                        scatter(x_range(i), 9, 50, 'r*');
                    end
                end
                if(save_file)
                    export_fig([OUTPATH '/0_FCD_lesion_feature_profile' ], '-m4', '-png'); close(gcf);
                end
            end
            
        end
        
        % 1.2 Group separation along the histological types: FCD IIa, IIb
        for separation_analysis_by_histo_types = 1
            
            FCD_IIa = strncmp(histo_type, 'FCDIIa', 6)';
            FCD_IIb = strncmp(histo_type, 'FCDIIb', 6)';
            
            sum(transmantle(FCD_IIa)==1)/length(find(FCD_IIa))
            sum(transmantle(FCD_IIb)==1)/length(find(FCD_IIb))
            
            p_type_IIa = [];
            for i = 1 : length(feature_set_idx)
                pat_data = mean(mean_z_lesion(feature_set_idx{i},FCD_IIa), 1);
                cont_data_temp = reshape(mean(mean_z_lesion_cont(feature_set_idx{i},:, :), 2), size(mean(mean_z_lesion_cont(feature_set_idx{i},:, :), 2), 1), size(mean(mean_z_lesion_cont(feature_set_idx{i},:, :), 2), 3));
                cont_data = mean(cont_data_temp, 1);
                
                [ p h stats ] = ranksum(pat_data, cont_data);
                
                p_type_IIa = [ p_type_IIa p ];
            end
            
            p_type_IIb = [];
            for i = 1 : length(feature_set_idx)
                pat_data = mean(mean_z_lesion(feature_set_idx{i}, FCD_IIb), 1);
                cont_data_temp = reshape(mean(mean_z_lesion_cont(feature_set_idx{i},:, :), 2), size(mean(mean_z_lesion_cont(feature_set_idx{i},:, :), 2), 1), size(mean(mean_z_lesion_cont(feature_set_idx{i},:, :), 2), 3));
                cont_data = mean(cont_data_temp, 1);
                
                [ p h stats ] = ranksum(pat_data, cont_data);
                p_type_IIb = [ p_type_IIb p ];
            end
            p_value_set_1 = [p_type_IIa p_type_IIb];
            
            p_type_IIa_vs_IIb_1 = [];
            include = ones(size(case_num_pat_temp, 1), 1);
            for i = 1 : length(feature_set_idx)
                [p_temp,h,stats] = ranksum(mean(mean_z_lesion(feature_set_idx{i}, FCD_IIa), 1), mean(mean_z_lesion(feature_set_idx{i}, FCD_IIb), 1));
                p_type_IIa_vs_IIb_1 = [ p_type_IIa_vs_IIb_1 p_temp ];
            end
            
            feature_profile_typeIIa = [
                ones(1, size(find(FCD_IIa), 2))*NaN
                mean(mean_z_lesion(2:5, FCD_IIa), 1)
                mean(mean_z_lesion([7:9 16], FCD_IIa), 1)
                mean(mean_z_lesion(12:15, FCD_IIa), 1)
                mean(mean_z_lesion(20:22, FCD_IIa), 1)
                mean(mean_z_lesion(17, FCD_IIa), 1)
                mean(mean_z_lesion(18, FCD_IIa), 1)
                mean(mean_z_lesion(19, FCD_IIa), 1)
                mean(mean_z_lesion(30:33, FCD_IIa), 1)
                mean(mean_z_lesion([35:37 44], FCD_IIa), 1)
                mean(mean_z_lesion(40:43, FCD_IIa), 1)
                mean(mean_z_lesion(45:47, FCD_IIa), 1)
                ]';
            feature_profile_typeIIa(feature_profile_typeIIa > 100) = 10;
            
            feature_profile_typeIIb = [
                ones(1, size(find(FCD_IIb), 2))*NaN
                mean(mean_z_lesion(2:5, FCD_IIb), 1)
                mean(mean_z_lesion([7:9 16], FCD_IIb), 1)
                mean(mean_z_lesion(12:15, FCD_IIb), 1)
                mean(mean_z_lesion(20:22, FCD_IIb), 1)
                mean(mean_z_lesion(17, FCD_IIb), 1)
                mean(mean_z_lesion(18, FCD_IIb), 1)
                mean(mean_z_lesion(19, FCD_IIb), 1)
                mean(mean_z_lesion(30:33, FCD_IIb), 1)
                mean(mean_z_lesion([35:37 44], FCD_IIb), 1)
                mean(mean_z_lesion(40:43, FCD_IIb), 1)
                mean(mean_z_lesion(45:47, FCD_IIb), 1)
                ]';
            feature_profile_typeIIb(feature_profile_typeIIb > 100) = 10;
            
            if(visualization)
                f = figure; set(gcf, 'Position', [587 500 1200 500]);
                hold on;
                linecolor = [0.8 0.8 0.8];
                hp1 = plot([0 52], [0 0]);
                hp2 = plot([0 52], [2 2], '--');
                hp3 = plot([0 52], [-2 -2], '--');
                set(hp1, 'Color', linecolor);
                set(hp2, 'Color', linecolor);
                set(hp3, 'Color', linecolor);
                set(gcf, 'Color', 'w');
                xtick = [0 3 6 9 12 16 19 22 27 30 33 36];
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
                
                xlim([2 max(xtick)+2]);
                set(gca, 'XTick', xtick+0.5);
                %% 0    3       6         9         12        16       19      22       27      30        33        36
                label_feature_profile_type = {'', 'Int', 'P-grad', 'T-grad', 'Sub Int', 'Thick', 'Curv', 'Depth', 'Int', 'P-grad', 'T-grad', 'Sub Int'};
                set(gca, 'XTickLabel', label_feature_profile_type);
                rotateticklabel(gca, 90);
                ylabel('z scores');
                ylim([-5 10]);
                set(gca, 'YTick', [-5 -2 0 2 5]);
                
                hold on;
                x_range = xtick(2:end);
                FDR_threshold = FDR_threshold_final_lesion_profiling;
                for i = 1 : size(p_type_IIa, 2)
                    if(p_type_IIa(i) <= FDR_threshold)
                        scatter(x_range(i), 9, 50, 'rv');
                    elseif(p_type_IIa(i) < 0.05 && p_type_IIa(i) > FDR_threshold)
                        scatter(x_range(i), 9, 50, 'r*');
                    end
                end
                
                for i = 1 : size(p_type_IIb, 2)
                    if(p_type_IIb(i) <= FDR_threshold)
                        scatter(x_range(i)+1, 9, 50, 'rv');
                    elseif(p_type_IIb(i) < 0.05 && p_type_IIb(i) > FDR_threshold)
                        scatter(x_range(i)+1, 9, 50, 'r*');
                    end
                end
                
                if(save_file)
                    export_fig([OUTPATH '/0_FCD_lesion_feature_profile_typeIIa_IIb' ], '-m4', '-png'); close(gcf);
                end
            end
            
        end
        
        % 2 FCD IIa vs. IIb - cortical depth analysis
        % 2-1. T1 intracortical intensity
        for depth_anlaysis_T1_cortical_intensity = 1
            
            excluded_case = { };
            included_cases = find(~ismember(case_num_pat_temp, excluded_case));
            FCD_IIa_temp = FCD_IIa' & included_cases;
            FCD_IIb_temp = FCD_IIb' & included_cases;
            
            feature_idx = [ 2:5 ];
            feature_profile_temp      = mean_z_lesion(feature_idx, included_cases);
            feature_profile_cont_temp = reshape(mean(mean_z_lesion_cont(feature_idx, :, :), 2), size(mean(mean_z_lesion_cont(feature_idx, :, :), 2), 1), size(mean(mean_z_lesion_cont(feature_idx, :, :), 2), 3));
            
            FCDIIa_feature = feature_profile_temp(:, FCD_IIa)';
            FCDIIb_feature = feature_profile_temp(:, FCD_IIb)';
            
            p_type_IIa = [];
            p_type_IIb = [];
            for j = 1 : length(feature_idx)
                
                pat_data = FCDIIa_feature(:, j);
                cont_data = feature_profile_cont_temp(j, :)';
                [ p h stats ] = ranksum(pat_data, cont_data);
                
                p_type_IIa = [ p_type_IIa p ];
                
                pat_data = FCDIIb_feature(:, j);
                cont_data = feature_profile_cont_temp(j, :)';
                [ p h stats ] = ranksum(pat_data, cont_data);
                
                p_type_IIb = [ p_type_IIb p ];
            end
            
            p_type_IIa_vs_IIb_2_1 = [];
            for j = 1 : length(feature_idx)
                [ p h stats ] = ranksum(FCDIIa_feature(:, j), FCDIIb_feature(:, j));
                p_type_IIa_vs_IIb_2_1 = [ p_type_IIa_vs_IIb_2_1 p ];
            end
            
            p_value_set_2_1 = [p_type_IIa p_type_IIb];
            
            if(visualization)
                x_lim = [0 3+(length(feature_idx)-1)*4];
                y_lim = [-4 4];
                
                visualization_graph_cortical_depth_profile(feature_idx, p_value_set_2_1, FCDIIa_feature, FCDIIb_feature, x_lim, y_lim, FDR_threshold_final_cortical_depth_profiling);
                
                if(save_file)
                    export_fig([OUTPATH '/1_cortical_depth_analysis_T1_intensity_typeII_a_vs_b_new' ], '-m3', '-png', '-painters'); close(gcf);
                end
            end
            
            p_value_set_2_1 = [p_value_set_2_1 p_type_IIa_vs_IIb_2_1];
            
            %% linear mixed effects model
            surf_idx_A = repmat((1:length(feature_idx))', sum(FCD_IIa_temp==1), 1);
            surf_idx_B = repmat((1:length(feature_idx))', sum(FCD_IIb_temp==1), 1);
            surf_idx = [ surf_idx_A; surf_idx_B ];
            group_A = repmat({'IIa'}, length(surf_idx_A), 1);
            group_B = repmat({'IIb'}, length(surf_idx_B), 1);
            group = [ group_A; group_B ];
            sub_A = repmat(case_num_pat(FCD_IIa_temp==1), 1,length(feature_idx))';
            sub_B = repmat(case_num_pat(FCD_IIb_temp==1), 1,length(feature_idx))';
            sub = [ sub_A(:); sub_B(:) ];
            Surf_idx = term(surf_idx);
            Group = term(group);
            Sub = term(sub);
            clear I
            
            %% 1) model of interests
            FCDIIa_feature_temp = FCDIIa_feature';
            FCDIIb_feature_temp = FCDIIb_feature';
            
            %% the effect of the surface level in whole patients
            Y = [ FCDIIa_feature_temp(:);  FCDIIb_feature_temp(:); ];
            model = 1 + Group + Surf_idx + Group*Surf_idx + random(sub) + I;
            slm = SurfStatLinMod(Y, model);
            slm = SurfStatT(slm, surf_idx);
            p = tcdf(slm.t, slm.df);
            p = min([p, 1-p]);
            p_value_set_2_1 = [ p_value_set_2_1 p];
            
            %% the effect of the surface level in patients with type IIa
            Y = [ FCDIIa_feature_temp(:); ];
            model = 1 + term(surf_idx_A) + random(sub_A(:)) + I;
            slm = SurfStatLinMod(Y, model);
            slm = SurfStatT(slm, surf_idx_A);
            p = tcdf(slm.t, slm.df);
            p = min([p, 1-p]);
            p_value_set_2_1 = [ p_value_set_2_1 p];
            
            %% the effect of the surface level in patients with type IIb
            Y = [ FCDIIb_feature_temp(:); ];
            model = 1 + term(surf_idx_B) + random(sub_B(:)) + I;
            slm = SurfStatLinMod(Y, model);
            slm = SurfStatT(slm, surf_idx_B);
            p = tcdf(slm.t, slm.df);
            p = min([p, 1-p]);
            p_value_set_2_1 = [ p_value_set_2_1 p];
            
            %% the interaction effect between surface level and group
            Y = [ FCDIIa_feature_temp(:);  FCDIIb_feature_temp(:); ];
            model = 1 + Group + Surf_idx + Group*Surf_idx + random(sub) + I;
            slm = SurfStatLinMod(Y, model);
            slm = SurfStatT(slm, surf_idx.*Group.IIa - surf_idx.*Group.IIb);
            p = tcdf(slm.t, slm.df);
            p = min([p, 1-p]);
            p_value_set_2_1 = [ p_value_set_2_1 p];
            
            if(visualization_linearmodel)
                CorticalDepthVisualization(feature_idx, FCD_IIa_temp, FCD_IIb_temp, FCDIIa_feature, FCDIIb_feature, [-5 5], 1, 0.05, p_type_IIa_vs_IIb_2_1/2, p_type_IIa, p_type_IIb, FDR_threshold_final_cortical_depth_profiling);
                if(save_file_linearmodel)
                    export_fig([OUTPATH '/1_cortical_depth_analysis_T1_intensity_typeII_a_vs_b_CI' ], '-m3', '-png', '-painters'); close(gcf);
                end
            end
            
        end
        
        % 2-2. T1 perpendicular gradient
        for depth_anlaysis_T1_perpendicular_gradient = 1
            
            excluded_case = { };
            included_cases = find(~ismember(case_num_pat_temp, excluded_case));
            FCD_IIa_temp = FCD_IIa' & included_cases;
            FCD_IIb_temp = FCD_IIb' & included_cases;
            
            feature_idx = [ 7:9 16 ] ;
            feature_profile_temp      = mean_z_lesion(feature_idx, included_cases);
            feature_profile_cont_temp = reshape(mean(mean_z_lesion_cont(feature_idx, :, :), 2), size(mean(mean_z_lesion_cont(feature_idx, :, :), 2), 1), size(mean(mean_z_lesion_cont(feature_idx, :, :), 2), 3));
            
            FCDIIa_feature = feature_profile_temp(:, FCD_IIa)';
            FCDIIb_feature = feature_profile_temp(:, FCD_IIb)';
            
            p_type_IIa = [];
            p_type_IIb = [];
            for j = 1 : length(feature_idx)
                
                pat_data = FCDIIa_feature(:, j);
                cont_data = feature_profile_cont_temp(j, :)';
                [ p h stats ] = ranksum(pat_data, cont_data);
                
                p_type_IIa = [ p_type_IIa 2*min(p, 1-p) ];
                
                pat_data = FCDIIb_feature(:, j);
                cont_data = feature_profile_cont_temp(j, :)';
                [ p h stats ] = ranksum(pat_data, cont_data);
                
                p_type_IIb = [ p_type_IIb 2*min(p, 1-p) ];
            end
            
            p_type_IIa_vs_IIb_2_2 = [];
            for j = 1 : length(feature_idx)
                [ p h stats ] = ranksum(FCDIIa_feature(:, j), FCDIIb_feature(:, j));
                p_type_IIa_vs_IIb_2_2 = [ p_type_IIa_vs_IIb_2_2 p ];
            end
            
            p_value_set_2_2 = [p_type_IIa p_type_IIb];
            
            if(visualization)
                x_lim = [0 3+(length(feature_idx)-1)*4];
                y_lim = [-6 4];
                
                visualization_graph_cortical_depth_profile(feature_idx, p_value_set_2_2, FCDIIa_feature, FCDIIb_feature, x_lim, y_lim, FDR_threshold_final_cortical_depth_profiling);
                
                if(save_file)
                    export_fig([OUTPATH '/2_cortical_depth_analysis_T1_p_gradient_typeII_a_vs_b_new' ], '-m3', '-png', '-painters'); close(gcf);
                end
            end
            
            p_value_set_2_2 = [p_value_set_2_2 p_type_IIa_vs_IIb_2_2];
            
            %% linear mixed effects model
            surf_idx_A = repmat((1:length(feature_idx))', sum(FCD_IIa_temp==1), 1);
            surf_idx_B = repmat((1:length(feature_idx))', sum(FCD_IIb_temp==1), 1);
            surf_idx = [ surf_idx_A; surf_idx_B ];
            group_A = repmat({'IIa'}, length(surf_idx_A), 1);
            group_B = repmat({'IIb'}, length(surf_idx_B), 1);
            group = [ group_A; group_B ];
            sub_A = repmat(case_num_pat(FCD_IIa_temp==1), 1,length(feature_idx))';
            sub_B = repmat(case_num_pat(FCD_IIb_temp==1), 1,length(feature_idx))';
            sub = [ sub_A(:); sub_B(:) ];
            Surf_idx = term(surf_idx);
            Group = term(group);
            Sub = term(sub);
            clear I
            
            %% 1) model of interests
            FCDIIa_feature_temp = FCDIIa_feature';
            FCDIIb_feature_temp = FCDIIb_feature';
            
            %% the effect of the surface level in whole patients
            Y = [ FCDIIa_feature_temp(:);  FCDIIb_feature_temp(:); ];
            model = 1 + Group + Surf_idx + Group*Surf_idx + random(sub) + I;
            slm = SurfStatLinMod(Y, model);
            slm = SurfStatT(slm, surf_idx);
            p = tcdf(slm.t, slm.df);
            p = min([p, 1-p])
            p_value_set_2_2 = [ p_value_set_2_2 p];
            
            %% the effect of the surface level in patients with type IIa
            Y = [ FCDIIa_feature_temp(:); ];
            model = 1 + term(surf_idx_A) + random(sub_A(:)) + I;
            slm = SurfStatLinMod(Y, model);
            slm = SurfStatT(slm, surf_idx_A);
            p = tcdf(slm.t, slm.df);
            p = min([p, 1-p])
            p_value_set_2_2 = [ p_value_set_2_2 p];
            
            %% the effect of the surface level in patients with type IIb
            Y = [ FCDIIb_feature_temp(:); ];
            model = 1 + term(surf_idx_B) + random(sub_B(:)) + I;
            slm = SurfStatLinMod(Y, model);
            slm = SurfStatT(slm, surf_idx_B);
            p = tcdf(slm.t, slm.df);
            p = min([p, 1-p])
            p_value_set_2_2 = [ p_value_set_2_2 p];
            
            %% the interaction effect between surface level and group
            Y = [ FCDIIa_feature_temp(:);  FCDIIb_feature_temp(:); ];
            model = 1 + Group + Surf_idx + Group*Surf_idx + random(sub) + I;
            slm = SurfStatLinMod(Y, model);
            slm = SurfStatT(slm, surf_idx.*Group.IIa - surf_idx.*Group.IIb);
            p = tcdf(slm.t, slm.df);
            p = min([p, 1-p])
            p_value_set_2_2 = [ p_value_set_2_2 p];
            
            if(visualization_linearmodel)
                CorticalDepthVisualization(feature_idx, FCD_IIa_temp, FCD_IIb_temp, FCDIIa_feature, FCDIIb_feature, [-5 6], 1, 0.05, p_type_IIa_vs_IIb_2_2/2, p_type_IIa, p_type_IIb, FDR_threshold_final_cortical_depth_profiling);
                if(save_file_linearmodel)
                    export_fig([OUTPATH '/2_cortical_depth_analysis_T1_p_gradient_typeII_a_vs_b_CI2' ], '-m3', '-png', '-painters'); close(gcf);
                end
            end
            
        end
        
        % 2-3. T1 tangential gradient
        for depth_anlaysis_T1_tangential_gradient = 1
            
            excluded_case = { };
            included_cases = find(~ismember(case_num_pat_temp, excluded_case));
            FCD_IIa_temp = FCD_IIa' & included_cases;
            FCD_IIb_temp = FCD_IIb' & included_cases;
            
            feature_idx = [ 12:15 ] ;
            feature_profile_temp      = mean_z_lesion(feature_idx, included_cases);
            feature_profile_cont_temp = reshape(mean(mean_z_lesion_cont(feature_idx, :, :), 2), size(mean(mean_z_lesion_cont(feature_idx, :, :), 2), 1), size(mean(mean_z_lesion_cont(feature_idx, :, :), 2), 3));
            
            FCDIIa_feature = feature_profile_temp(:, FCD_IIa)';
            FCDIIb_feature = feature_profile_temp(:, FCD_IIb)';
            
            p_type_IIa = [];
            p_type_IIb = [];
            for j = 1 : length(feature_idx)
                
                pat_data = FCDIIa_feature(:, j);
                cont_data = feature_profile_cont_temp(j, :)';
                [ p h stats ] = ranksum(pat_data, cont_data);
                
                p_type_IIa = [ p_type_IIa 2*min(p, 1-p) ];
                
                pat_data = FCDIIb_feature(:, j);
                cont_data = feature_profile_cont_temp(j, :)';
                [ p h stats ] = ranksum(pat_data, cont_data);
                
                p_type_IIb = [ p_type_IIb 2*min(p, 1-p) ];
            end
            
            p_type_IIa_vs_IIb_2_3 = [];
            for j = 1 : length(feature_idx)
                [ p h stats ] = ranksum(FCDIIa_feature(:, j), FCDIIb_feature(:, j));
                p_type_IIa_vs_IIb_2_3 = [ p_type_IIa_vs_IIb_2_3 p ];
            end
            p_value_set_2_3 = [p_type_IIa p_type_IIb];
            
            if(visualization)
                x_lim = [0 3+(length(feature_idx)-1)*4];
                y_lim = [-4 3];
                
                visualization_graph_cortical_depth_profile(feature_idx, p_value_set_2_3, FCDIIa_feature, FCDIIb_feature, x_lim, y_lim, FDR_threshold_final_cortical_depth_profiling);
                
                if(save_file)
                    export_fig([OUTPATH '/3_cortical_depth_analysis_T1_t_gradient_typeII_a_vs_b_new' ], '-m3', '-png', '-painters'); close(gcf);
                end
            end
            
            p_value_set_2_3 = [p_value_set_2_3 p_type_IIa_vs_IIb_2_3];
            
            %% linear mixed effects model
            surf_idx_A = repmat((1:length(feature_idx))', sum(FCD_IIa_temp==1), 1);
            surf_idx_B = repmat((1:length(feature_idx))', sum(FCD_IIb_temp==1), 1);
            surf_idx = [ surf_idx_A; surf_idx_B ];
            group_A = repmat({'IIa'}, length(surf_idx_A), 1);
            group_B = repmat({'IIb'}, length(surf_idx_B), 1);
            group = [ group_A; group_B ];
            sub_A = repmat(case_num_pat(FCD_IIa_temp==1), 1,length(feature_idx))';
            sub_B = repmat(case_num_pat(FCD_IIb_temp==1), 1,length(feature_idx))';
            sub = [ sub_A(:); sub_B(:) ];
            Surf_idx = term(surf_idx);
            Group = term(group);
            Sub = term(sub);
            clear I
            
            %% 1) model of interests
            FCDIIa_feature_temp = FCDIIa_feature';
            FCDIIb_feature_temp = FCDIIb_feature';
            
            %% the effect of the surface level in whole patients
            Y = [ FCDIIa_feature_temp(:);  FCDIIb_feature_temp(:); ];
            model = 1 + Group + Surf_idx + Group*Surf_idx + random(sub) + I;
            slm = SurfStatLinMod(Y, model);
            slm = SurfStatT(slm, surf_idx);
            p = tcdf(slm.t, slm.df);
            p = min([p, 1-p])
            p_value_set_2_3 = [ p_value_set_2_3 p];
            
            %% the effect of the surface level in patients with type IIa
            Y = [ FCDIIa_feature_temp(:); ];
            model = 1 + term(surf_idx_A) + random(sub_A(:)) + I;
            slm = SurfStatLinMod(Y, model);
            slm = SurfStatT(slm, surf_idx_A);
            p = tcdf(slm.t, slm.df);
            p = min([p, 1-p])
            p_value_set_2_3 = [ p_value_set_2_3 p];
            
            %% the effect of the surface level in patients with type IIb
            Y = [ FCDIIb_feature_temp(:); ];
            model = 1 + term(surf_idx_B) + random(sub_B(:)) + I;
            slm = SurfStatLinMod(Y, model);
            slm = SurfStatT(slm, surf_idx_B);
            p = tcdf(slm.t, slm.df);
            p = min([p, 1-p])
            p_value_set_2_3 = [ p_value_set_2_3 p];
            
            %% the interaction effect between surface level and group
            Y = [ FCDIIa_feature_temp(:);  FCDIIb_feature_temp(:); ];
            model = 1 + Group + Surf_idx + Group*Surf_idx + random(sub) + I;
            slm = SurfStatLinMod(Y, model);
            slm = SurfStatT(slm, surf_idx.*Group.IIa - surf_idx.*Group.IIb);
            p = tcdf(slm.t, slm.df);
            p = min([p, 1-p])
            p_value_set_2_3 = [ p_value_set_2_3 p];
            
            if(visualization_linearmodel)
                CorticalDepthVisualization(feature_idx, FCD_IIa_temp, FCD_IIb_temp, FCDIIa_feature, FCDIIb_feature, [-5 5], 1, 0.05, p_type_IIa_vs_IIb_2_3/2, p_type_IIa, p_type_IIb, FDR_threshold_final_cortical_depth_profiling);
                if(save_file_linearmodel)
                    export_fig([OUTPATH '/3_cortical_depth_analysis_T1_t_gradient_typeII_a_vs_b_CI2' ], '-m3', '-png', '-painters'); close(gcf);
                end
            end
            
        end
        
        % 2-4. T1 subcortical intensity
        for depth_anlaysis_T1_subcortical_intensity = 1
            
            excluded_case = { };
            included_cases = find(~ismember(case_num_pat_temp, excluded_case));
            FCD_IIa_temp = FCD_IIa' & included_cases;
            FCD_IIb_temp = FCD_IIb' & included_cases;
            
            feature_idx = [ 20:22 ] ;
            feature_profile_temp      = mean_z_lesion(feature_idx, included_cases);
            feature_profile_cont_temp = reshape(mean(mean_z_lesion_cont(feature_idx, :, :), 2), size(mean(mean_z_lesion_cont(feature_idx, :, :), 2), 1), size(mean(mean_z_lesion_cont(feature_idx, :, :), 2), 3));
            
            FCDIIa_feature = feature_profile_temp(:, FCD_IIa)';
            FCDIIb_feature = feature_profile_temp(:, FCD_IIb)';
            
            p_type_IIa = [];
            p_type_IIb = [];
            for j = 1 : length(feature_idx)
                
                pat_data = FCDIIa_feature(:, j);
                cont_data = feature_profile_cont_temp(j, :)';
                [ p h stats ] = ranksum(pat_data, cont_data);
                
                p_type_IIa = [ p_type_IIa 2*min(p, 1-p) ];
                
                pat_data = FCDIIb_feature(:, j);
                cont_data = feature_profile_cont_temp(j, :)';
                [ p h stats ] = ranksum(pat_data, cont_data);
                
                p_type_IIb = [ p_type_IIb 2*min(p, 1-p) ];
            end
            
            p_type_IIa_vs_IIb_2_4 = [];
            for j = 1 : length(feature_idx)
                [ p h stats ] = ranksum(FCDIIa_feature(:, j), FCDIIb_feature(:, j));
                p_type_IIa_vs_IIb_2_4 = [ p_type_IIa_vs_IIb_2_4 p ];
            end
            
            p_value_set_2_4 = [p_type_IIa p_type_IIb];
            
            if(visualization)
                x_lim = [0 3+(length(feature_idx)-1)*4];
                y_lim = [-6 3];
                
                visualization_graph_cortical_depth_profile(feature_idx, p_value_set_2_4, FCDIIa_feature, FCDIIb_feature, x_lim, y_lim, FDR_threshold_final_cortical_depth_profiling);
                if(save_file)
                    export_fig([OUTPATH '/4_cortical_depth_analysis_T1_subintensity_typeII_a_vs_b_new' ], '-m3', '-png', '-painters'); close(gcf);
                end
            end
            
            p_value_set_2_4 = [p_value_set_2_4 p_type_IIa_vs_IIb_2_4];
            
            %% linear mixed effects model
            surf_idx_A = repmat((1:length(feature_idx))', sum(FCD_IIa_temp==1), 1);
            surf_idx_B = repmat((1:length(feature_idx))', sum(FCD_IIb_temp==1), 1);
            surf_idx = [ surf_idx_A; surf_idx_B ];
            group_A = repmat({'IIa'}, length(surf_idx_A), 1);
            group_B = repmat({'IIb'}, length(surf_idx_B), 1);
            group = [ group_A; group_B ];
            sub_A = repmat(case_num_pat(FCD_IIa_temp==1), 1,length(feature_idx))';
            sub_B = repmat(case_num_pat(FCD_IIb_temp==1), 1,length(feature_idx))';
            sub = [ sub_A(:); sub_B(:) ];
            Surf_idx = term(surf_idx);
            Group = term(group);
            Sub = term(sub);
            clear I
            
            %% 1) model of interests
            FCDIIa_feature_temp = FCDIIa_feature';
            FCDIIb_feature_temp = FCDIIb_feature';
            
            %% the effect of the surface level in whole patients
            Y = [ FCDIIa_feature_temp(:);  FCDIIb_feature_temp(:); ];
            model = 1 + Group + Surf_idx + Group*Surf_idx + random(sub) + I;
            slm = SurfStatLinMod(Y, model);
            slm = SurfStatT(slm, surf_idx);
            p = tcdf(slm.t, slm.df);
            p = min([p, 1-p])
            p_value_set_2_4 = [ p_value_set_2_4 p];
            
            %% the effect of the surface level in patients with type IIa
            Y = [ FCDIIa_feature_temp(:); ];
            model = 1 + term(surf_idx_A) + random(sub_A(:)) + I;
            slm = SurfStatLinMod(Y, model);
            slm = SurfStatT(slm, surf_idx_A);
            p = tcdf(slm.t, slm.df);
            p = min([p, 1-p])
            p_value_set_2_4 = [ p_value_set_2_4 p];
            
            %% the effect of the surface level in patients with type IIb
            Y = [ FCDIIb_feature_temp(:); ];
            model = 1 + term(surf_idx_B) + random(sub_B(:)) + I;
            slm = SurfStatLinMod(Y, model);
            slm = SurfStatT(slm, surf_idx_B);
            p = tcdf(slm.t, slm.df);
            p = min([p, 1-p])
            p_value_set_2_4 = [ p_value_set_2_4 p];
            
            %% the interaction effect between surface level and group
            Y = [ FCDIIa_feature_temp(:);  FCDIIb_feature_temp(:); ];
            model = 1 + Group + Surf_idx + Group*Surf_idx + random(sub) + I;
            slm = SurfStatLinMod(Y, model);
            slm = SurfStatT(slm, surf_idx.*Group.IIa - surf_idx.*Group.IIb);
            p = tcdf(slm.t, slm.df);
            p = min([p, 1-p])
            p_value_set_2_4 = [ p_value_set_2_4 p];
            
            if(visualization_linearmodel)
                CorticalDepthVisualization(feature_idx, FCD_IIa_temp, FCD_IIb_temp, FCDIIa_feature, FCDIIb_feature, [-6 6], 1, 0.05, p_type_IIa_vs_IIb_2_4/2, p_type_IIa, p_type_IIb, FDR_threshold_final_cortical_depth_profiling);
                if(save_file_linearmodel)
                    export_fig([OUTPATH '/4_cortical_depth_analysis_T1_subintensity_typeII_a_vs_b_CI2' ], '-m3', '-png', '-painters'); close(gcf);
                end
            end
            
        end
        
        % 2-5. FLAIR intracortical intensity
        for depth_anlaysis_FLAIR_cortical_intensity = 1
            
            excluded_case = { };
            included_cases = find(~ismember(case_num_pat_temp, excluded_case));
            included_cases = ~ismember(case_num_pat_temp, excluded_case);
            FCD_IIa_temp = FCD_IIa' & included_cases;
            FCD_IIb_temp = FCD_IIb' & included_cases;
            
            feature_idx = [ 30:33 ] ;
            feature_profile_temp      = mean_z_lesion(feature_idx, included_cases);
            feature_profile_cont_temp = reshape(mean(mean_z_lesion_cont(feature_idx, :, :), 2), size(mean(mean_z_lesion_cont(feature_idx, :, :), 2), 1), size(mean(mean_z_lesion_cont(feature_idx, :, :), 2), 3));
            
            FCDIIa_feature = feature_profile_temp(:, FCD_IIa)';
            FCDIIb_feature = feature_profile_temp(:, FCD_IIb)';
            
            p_type_IIa = [];
            p_type_IIb = [];
            for j = 1 : length(feature_idx)
                
                pat_data = FCDIIa_feature(:, j);
                cont_data = feature_profile_cont_temp(j, :)';
                [ p h stats ] = ranksum(pat_data, cont_data);
                
                p_type_IIa = [ p_type_IIa 2*min(p, 1-p) ];
                
                pat_data = FCDIIb_feature(:, j);
                cont_data = feature_profile_cont_temp(j, :)';
                [ p h stats ] = ranksum(pat_data, cont_data);
                
                p_type_IIb = [ p_type_IIb 2*min(p, 1-p) ];
            end
            
            p_type_IIa_vs_IIb_2_7 = [];
            for j = 1 : length(feature_idx)
                [ p h stats ] = ranksum(FCDIIa_feature(:, j), FCDIIb_feature(:, j));
                p_type_IIa_vs_IIb_2_7 = [ p_type_IIa_vs_IIb_2_7 p ];
            end
            
            p_value_set_2_7 = [p_type_IIa p_type_IIb];
            
            if(visualization)
                x_lim = [0 3+(length(feature_idx)-1)*4];
                y_lim = [-3 6];
                
                visualization_graph_cortical_depth_profile(feature_idx, p_value_set_2_7, FCDIIa_feature, FCDIIb_feature, x_lim, y_lim, FDR_threshold_final_cortical_depth_profiling);
                
                if(save_file)
                    export_fig([OUTPATH '/5_cortical_depth_analysis_FLAIR_intensity_typeII_a_vs_b_new' ], '-m3', '-png', '-painters'); close(gcf);
                end
            end
            
            p_value_set_2_7 = [p_value_set_2_7 p_type_IIa_vs_IIb_2_7];
            
            %% linear mixed effects model
            surf_idx_A = repmat((1:length(feature_idx))', sum(FCD_IIa_temp==1), 1);
            surf_idx_B = repmat((1:length(feature_idx))', sum(FCD_IIb_temp==1), 1);
            surf_idx = [ surf_idx_A; surf_idx_B ];
            group_A = repmat({'IIa'}, length(surf_idx_A), 1);
            group_B = repmat({'IIb'}, length(surf_idx_B), 1);
            group = [ group_A; group_B ];
            sub_A = repmat(case_num_pat(FCD_IIa_temp==1), 1,length(feature_idx))';
            sub_B = repmat(case_num_pat(FCD_IIb_temp==1), 1,length(feature_idx))';
            sub = [ sub_A(:); sub_B(:) ];
            Surf_idx = term(surf_idx);
            Group = term(group);
            Sub = term(sub);
            clear I
            
            %% 1) model of interests
            FCDIIa_feature_temp = FCDIIa_feature';
            FCDIIb_feature_temp = FCDIIb_feature';
            
            %% the effect of the surface level in whole patients
            Y = [ FCDIIa_feature_temp(:);  FCDIIb_feature_temp(:); ];
            model = 1 + Group + Surf_idx + Group*Surf_idx + random(sub) + I;
            slm = SurfStatLinMod(Y, model);
            slm = SurfStatT(slm, surf_idx);
            p = tcdf(slm.t, slm.df);
            p = min([p, 1-p])
            p_value_set_2_7 = [ p_value_set_2_7 p];
            
            %% the effect of the surface level in patients with type IIa
            Y = [ FCDIIa_feature_temp(:); ];
            model = 1 + term(surf_idx_A) + random(sub_A(:)) + I;
            slm = SurfStatLinMod(Y, model);
            slm = SurfStatT(slm, surf_idx_A);
            p = tcdf(slm.t, slm.df);
            p = min([p, 1-p])
            p_value_set_2_7 = [ p_value_set_2_7 p];
            
            %% the effect of the surface level in patients with type IIb
            Y = [ FCDIIb_feature_temp(:); ];
            model = 1 + term(surf_idx_B) + random(sub_B(:)) + I;
            slm = SurfStatLinMod(Y, model);
            slm = SurfStatT(slm, surf_idx_B);
            p = tcdf(slm.t, slm.df);
            p = min([p, 1-p])
            p_value_set_2_7 = [ p_value_set_2_7 p];
            
            %% the interaction effect between surface level and group
            Y = [ FCDIIa_feature_temp(:);  FCDIIb_feature_temp(:); ];
            model = 1 + Group + Surf_idx + Group*Surf_idx + random(sub) + I;
            slm = SurfStatLinMod(Y, model);
            slm = SurfStatT(slm, surf_idx.*Group.IIa - surf_idx.*Group.IIb);
            p = tcdf(slm.t, slm.df);
            p = min([p, 1-p])
            p_value_set_2_7 = [ p_value_set_2_7 p];
            
            if(visualization_linearmodel)
                CorticalDepthVisualization(feature_idx, FCD_IIa_temp, FCD_IIb_temp, FCDIIa_feature, FCDIIb_feature, [-3 7], 1, 0.05, p_type_IIa_vs_IIb_2_7/2, p_type_IIa, p_type_IIb, FDR_threshold_final_cortical_depth_profiling);
                if(save_file_linearmodel)
                    export_fig([OUTPATH '/5_cortical_depth_analysis_FLAIR_intensity_typeII_a_vs_b_CI2' ], '-m3', '-png', '-painters'); close(gcf);
                end
            end
            
        end
        
        % 2-6. FLAIR perpendicular gradient
        for depth_anlaysis_FLAIR_perpendicular_gradient = 1
            
            excluded_case = { };
            included_cases = find(~ismember(case_num_pat_temp, excluded_case));
            FCD_IIa_temp = FCD_IIa' & included_cases;
            FCD_IIb_temp = FCD_IIb' & included_cases;
            
            feature_idx = [ 35:37 44] ;
            feature_profile_temp      = mean_z_lesion(feature_idx, included_cases);
            feature_profile_cont_temp = reshape(mean(mean_z_lesion_cont(feature_idx, :, :), 2), size(mean(mean_z_lesion_cont(feature_idx, :, :), 2), 1), size(mean(mean_z_lesion_cont(feature_idx, :, :), 2), 3));
            
            FCDIIa_feature = feature_profile_temp(:, FCD_IIa)';
            FCDIIb_feature = feature_profile_temp(:, FCD_IIb)';
            
            p_type_IIa = [];
            p_type_IIb = [];
            for j = 1 : length(feature_idx)
                
                pat_data = FCDIIa_feature(:, j);
                cont_data = feature_profile_cont_temp(j, :)';
                [ p h stats ] = ranksum(pat_data, cont_data);
                
                p_type_IIa = [ p_type_IIa 2*min(p, 1-p) ];
                
                pat_data = FCDIIb_feature(:, j);
                cont_data = feature_profile_cont_temp(j, :)';
                [ p h stats ] = ranksum(pat_data, cont_data);
                
                p_type_IIb = [ p_type_IIb 2*min(p, 1-p) ];
            end
            
            p_type_IIa_vs_IIb_2_8 = [];
            for j = 1 : length(feature_idx)
                [ p h stats ] = ranksum(FCDIIa_feature(:, j), FCDIIb_feature(:, j));
                p_type_IIa_vs_IIb_2_8 = [ p_type_IIa_vs_IIb_2_8 p ];
            end
            p_value_set_2_8 = [p_type_IIa p_type_IIb];
            
            if(visualization)
                x_lim = [0 3+(length(feature_idx)-1)*4];
                y_lim = [-4 3];
                
                visualization_graph_cortical_depth_profile(feature_idx, p_value_set_2_8, FCDIIa_feature, FCDIIb_feature, x_lim, y_lim, FDR_threshold_final_cortical_depth_profiling);
                
                if(save_file)
                    export_fig([OUTPATH '/6_cortical_depth_analysis_FLAIR_p_gradient_typeII_a_vs_b_new' ], '-m3', '-png', '-painters'); close(gcf);
                end
            end
            
            p_value_set_2_8 = [p_value_set_2_8 p_type_IIa_vs_IIb_2_8];
            
            %% linear mixed effects model
            surf_idx_A = repmat((1:length(feature_idx))', sum(FCD_IIa_temp==1), 1);
            surf_idx_B = repmat((1:length(feature_idx))', sum(FCD_IIb_temp==1), 1);
            surf_idx = [ surf_idx_A; surf_idx_B ];
            group_A = repmat({'IIa'}, length(surf_idx_A), 1);
            group_B = repmat({'IIb'}, length(surf_idx_B), 1);
            group = [ group_A; group_B ];
            sub_A = repmat(case_num_pat(FCD_IIa_temp==1), 1,length(feature_idx))';
            sub_B = repmat(case_num_pat(FCD_IIb_temp==1), 1,length(feature_idx))';
            sub = [ sub_A(:); sub_B(:) ];
            Surf_idx = term(surf_idx);
            Group = term(group);
            Sub = term(sub);
            clear I
            
            %% 1) model of interests
            FCDIIa_feature_temp = FCDIIa_feature';
            FCDIIb_feature_temp = FCDIIb_feature';
            
            %% the effect of the surface level in whole patients
            Y = [ FCDIIa_feature_temp(:);  FCDIIb_feature_temp(:); ];
            model = 1 + Group + Surf_idx + Group*Surf_idx + random(sub) + I;
            slm = SurfStatLinMod(Y, model);
            slm = SurfStatT(slm, surf_idx);
            p = tcdf(slm.t, slm.df);
            p = min([p, 1-p])
            p_value_set_2_8 = [ p_value_set_2_8 p];
            
            %% the effect of the surface level in patients with type IIa
            Y = [ FCDIIa_feature_temp(:); ];
            model = 1 + term(surf_idx_A) + random(sub_A(:)) + I;
            slm = SurfStatLinMod(Y, model);
            slm = SurfStatT(slm, surf_idx_A);
            p = tcdf(slm.t, slm.df);
            p = min([p, 1-p])
            p_value_set_2_8 = [ p_value_set_2_8 p];
            
            %% the effect of the surface level in patients with type IIb
            Y = [ FCDIIb_feature_temp(:); ];
            model = 1 + term(surf_idx_B) + random(sub_B(:)) + I;
            slm = SurfStatLinMod(Y, model);
            slm = SurfStatT(slm, surf_idx_B);
            p = tcdf(slm.t, slm.df);
            p = min([p, 1-p])
            p_value_set_2_8 = [ p_value_set_2_8 p];
            
            %% the interaction effect between surface level and group
            Y = [ FCDIIa_feature_temp(:);  FCDIIb_feature_temp(:); ];
            model = 1 + Group + Surf_idx + Group*Surf_idx + random(sub) + I;
            slm = SurfStatLinMod(Y, model);
            slm = SurfStatT(slm, surf_idx.*Group.IIa - surf_idx.*Group.IIb);
            p = tcdf(slm.t, slm.df);
            p = min([p, 1-p])
            p_value_set_2_8 = [ p_value_set_2_8 p];
            
            if(visualization_linearmodel)
                CorticalDepthVisualization(feature_idx, FCD_IIa_temp, FCD_IIb_temp, FCDIIa_feature, FCDIIb_feature, [-5 5], 1, 0.05, p_type_IIa_vs_IIb_2_8/2, p_type_IIa, p_type_IIb, FDR_threshold_final_cortical_depth_profiling);
                if(save_file_linearmodel)
                    export_fig([OUTPATH '/6_cortical_depth_analysis_FLAIR_p_gradient_typeII_a_vs_b_CI' ], '-m3', '-png', '-painters'); close(gcf);
                end
            end
            
        end
        
        % 2-7. FLAIR tangential gradient
        for depth_anlaysis_FLAIR_tangential_gradient = 1
            
            excluded_case = { };
            included_cases = find(~ismember(case_num_pat_temp, excluded_case));
            
            feature_idx = [ 40:43 ] ;
            feature_profile_temp      = mean_z_lesion(feature_idx, included_cases);
            feature_profile_cont_temp = reshape(mean(mean_z_lesion_cont(feature_idx, :, :), 2), size(mean(mean_z_lesion_cont(feature_idx, :, :), 2), 1), size(mean(mean_z_lesion_cont(feature_idx, :, :), 2), 3));
            
            FCDIIa_feature = feature_profile_temp(:, FCD_IIa)';
            FCDIIb_feature = feature_profile_temp(:, FCD_IIb)';
            
            p_type_IIa = [];
            p_type_IIb = [];
            for j = 1 : length(feature_idx)
                
                pat_data = FCDIIa_feature(:, j);
                cont_data = feature_profile_cont_temp(j, :)';
                [ p h stats ] = ranksum(pat_data, cont_data);
                
                p_type_IIa = [ p_type_IIa 2*min(p, 1-p) ];
                
                pat_data = FCDIIb_feature(:, j);
                cont_data = feature_profile_cont_temp(j, :)';
                [ p h stats ] = ranksum(pat_data, cont_data);
                
                p_type_IIb = [ p_type_IIb 2*min(p, 1-p) ];
            end
            
            p_type_IIa_vs_IIb_2_9 = [];
            for j = 1 : length(feature_idx)
                [ p h stats ] = ranksum(FCDIIa_feature(:, j), FCDIIb_feature(:, j));
                p_type_IIa_vs_IIb_2_9 = [ p_type_IIa_vs_IIb_2_9 p ];
            end
            
            p_value_set_2_9 = [p_type_IIa p_type_IIb];
            
            if(visualization)
                x_lim = [0 3+(length(feature_idx)-1)*4];
                y_lim = [-3 4];
                
                visualization_graph_cortical_depth_profile(feature_idx, p_value_set_2_9, FCDIIa_feature, FCDIIb_feature, x_lim, y_lim, FDR_threshold_final_cortical_depth_profiling);
                
                if(save_file)
                    export_fig([OUTPATH '/7_cortical_depth_analysis_FLAIR_t_gradient_typeII_a_vs_b_new' ], '-m3', '-png', '-painters'); close(gcf);
                end
            end
            
            p_value_set_2_9 = [p_value_set_2_9 p_type_IIa_vs_IIb_2_9];
            
            %% linear mixed effects model
            surf_idx_A = repmat((1:length(feature_idx))', sum(FCD_IIa_temp==1), 1);
            surf_idx_B = repmat((1:length(feature_idx))', sum(FCD_IIb_temp==1), 1);
            surf_idx = [ surf_idx_A; surf_idx_B ];
            group_A = repmat({'IIa'}, length(surf_idx_A), 1);
            group_B = repmat({'IIb'}, length(surf_idx_B), 1);
            group = [ group_A; group_B ];
            sub_A = repmat(case_num_pat(FCD_IIa_temp==1), 1,length(feature_idx))';
            sub_B = repmat(case_num_pat(FCD_IIb_temp==1), 1,length(feature_idx))';
            sub = [ sub_A(:); sub_B(:) ];
            Surf_idx = term(surf_idx);
            Group = term(group);
            Sub = term(sub);
            clear I
            
            %% 1) model of interests
            FCDIIa_feature_temp = FCDIIa_feature';
            FCDIIb_feature_temp = FCDIIb_feature';
            
            %% the effect of the surface level in whole patients
            Y = [ FCDIIa_feature_temp(:);  FCDIIb_feature_temp(:); ];
            model = 1 + Group + Surf_idx + Group*Surf_idx + random(sub) + I;
            slm = SurfStatLinMod(Y, model);
            slm = SurfStatT(slm, surf_idx);
            p = tcdf(slm.t, slm.df);
            p = min([p, 1-p])
            p_value_set_2_9 = [ p_value_set_2_9 p];
            
            %% the effect of the surface level in patients with type IIa
            Y = [ FCDIIa_feature_temp(:); ];
            model = 1 + term(surf_idx_A) + random(sub_A(:)) + I;
            slm = SurfStatLinMod(Y, model);
            slm = SurfStatT(slm, surf_idx_A);
            p = tcdf(slm.t, slm.df);
            p = min([p, 1-p])
            p_value_set_2_9 = [ p_value_set_2_9 p];
            
            %% the effect of the surface level in patients with type IIb
            Y = [ FCDIIb_feature_temp(:); ];
            model = 1 + term(surf_idx_B) + random(sub_B(:)) + I;
            slm = SurfStatLinMod(Y, model);
            slm = SurfStatT(slm, surf_idx_B);
            p = tcdf(slm.t, slm.df);
            p = min([p, 1-p])
            p_value_set_2_9 = [ p_value_set_2_9 p];
            
            %% the interaction effect between surface level and group
            Y = [ FCDIIa_feature_temp(:);  FCDIIb_feature_temp(:); ];
            model = 1 + Group + Surf_idx + Group*Surf_idx + random(sub) + I;
            slm = SurfStatLinMod(Y, model);
            slm = SurfStatT(slm, surf_idx.*Group.IIa - surf_idx.*Group.IIb);
            p = tcdf(slm.t, slm.df);
            p = min([p, 1-p])
            p_value_set_2_9 = [ p_value_set_2_9 p];
            
            if(visualization_linearmodel)
                CorticalDepthVisualization(feature_idx, FCD_IIa_temp, FCD_IIb_temp, FCDIIa_feature, FCDIIb_feature, [-3 3], 1, 0.05, p_type_IIa_vs_IIb_2_9/2, p_type_IIa, p_type_IIb, FDR_threshold_final_cortical_depth_profiling);
                if(save_file_linearmodel)
                    export_fig([OUTPATH '/7_cortical_depth_analysis_FLAIR_t_gradient_typeII_a_vs_b_CI' ], '-m3', '-png', '-painters'); close(gcf);
                end
            end
            
        end
        
        % 2-8. FLAIR subcortical intensity
        for depth_anlaysis_FLAIR_subcortical_intensity = 1
            
            excluded_case = { };
            included_cases = find(~ismember(case_num_pat_temp, excluded_case));
            FCD_IIa_temp = FCD_IIa' & included_cases;
            FCD_IIb_temp = FCD_IIb' & included_cases;
            
            feature_idx = [ 45:47 ] ;
            feature_profile_temp      = mean_z_lesion(feature_idx, included_cases);
            feature_profile_cont_temp = reshape(mean(mean_z_lesion_cont(feature_idx, :, :), 2), size(mean(mean_z_lesion_cont(feature_idx, :, :), 2), 1), size(mean(mean_z_lesion_cont(feature_idx, :, :), 2), 3));
            
            FCDIIa_feature = feature_profile_temp(:, FCD_IIa)';
            FCDIIb_feature = feature_profile_temp(:, FCD_IIb)';
            
            p_type_IIa = [];
            p_type_IIb = [];
            for j = 1 : length(feature_idx)
                
                pat_data = FCDIIa_feature(:, j);
                cont_data = feature_profile_cont_temp(j, :)';
                [ p h stats ] = ranksum(pat_data, cont_data);
                
                p_type_IIa = [ p_type_IIa 2*min(p, 1-p) ];
                
                pat_data = FCDIIb_feature(:, j);
                cont_data = feature_profile_cont_temp(j, :)';
                [ p h stats ] = ranksum(pat_data, cont_data);
                
                p_type_IIb = [ p_type_IIb 2*min(p, 1-p) ];
            end
            
            p_type_IIa_vs_IIb_2_10 = [];
            for j = 1 : length(feature_idx)
                [ p h stats ] = ranksum(FCDIIa_feature(:, j), FCDIIb_feature(:, j));
                p_type_IIa_vs_IIb_2_10 = [ p_type_IIa_vs_IIb_2_10 p ];
            end
            p_value_set_2_10 = [p_type_IIa p_type_IIb];
            
            if(visualization)
                x_lim = [0 3+(length(feature_idx)-1)*4];
                y_lim = [-2 6];
                
                visualization_graph_cortical_depth_profile(feature_idx, p_value_set_2_10, FCDIIa_feature, FCDIIb_feature, x_lim, y_lim, FDR_threshold_final_cortical_depth_profiling);
                
                if(save_file)
                    export_fig([OUTPATH '/8_cortical_depth_analysis_FLAIR_sub_intensity_typeII_a_vs_b_new' ], '-m3', '-png', '-painters'); close(gcf);
                end
            end
            
            p_value_set_2_10 = [p_value_set_2_10 p_type_IIa_vs_IIb_2_10];
            
            %% linear mixed effects model
            surf_idx_A = repmat((1:length(feature_idx))', sum(FCD_IIa_temp==1), 1);
            surf_idx_B = repmat((1:length(feature_idx))', sum(FCD_IIb_temp==1), 1);
            surf_idx = [ surf_idx_A; surf_idx_B ];
            group_A = repmat({'IIa'}, length(surf_idx_A), 1);
            group_B = repmat({'IIb'}, length(surf_idx_B), 1);
            group = [ group_A; group_B ];
            sub_A = repmat(case_num_pat(FCD_IIa_temp==1), 1,length(feature_idx))';
            sub_B = repmat(case_num_pat(FCD_IIb_temp==1), 1,length(feature_idx))';
            sub = [ sub_A(:); sub_B(:) ];
            Surf_idx = term(surf_idx);
            Group = term(group);
            Sub = term(sub);
            clear I
            
            %% 1) model of interests
            FCDIIa_feature_temp = FCDIIa_feature';
            FCDIIb_feature_temp = FCDIIb_feature';
            
            %% the effect of the surface level in whole patients
            Y = [ FCDIIa_feature_temp(:);  FCDIIb_feature_temp(:); ];
            model = 1 + Group + Surf_idx + Group*Surf_idx + random(sub) + I;
            slm = SurfStatLinMod(Y, model);
            slm = SurfStatT(slm, surf_idx);
            p = tcdf(slm.t, slm.df);
            p = min([p, 1-p])
            p_value_set_2_10 = [ p_value_set_2_10 p];
            
            %% the effect of the surface level in patients with type IIa
            Y = [ FCDIIa_feature_temp(:); ];
            model = 1 + term(surf_idx_A) + random(sub_A(:)) + I;
            slm = SurfStatLinMod(Y, model);
            slm = SurfStatT(slm, surf_idx_A);
            p = tcdf(slm.t, slm.df);
            p = min([p, 1-p])
            p_value_set_2_10 = [ p_value_set_2_10 p];
            
            %% the effect of the surface level in patients with type IIb
            Y = [ FCDIIb_feature_temp(:); ];
            model = 1 + term(surf_idx_B) + random(sub_B(:)) + I;
            slm = SurfStatLinMod(Y, model);
            slm = SurfStatT(slm, surf_idx_B);
            p = tcdf(slm.t, slm.df);
            p = min([p, 1-p])
            p_value_set_2_10 = [ p_value_set_2_10 p];
            
            %% the interaction effect between surface level and group
            Y = [ FCDIIa_feature_temp(:);  FCDIIb_feature_temp(:); ];
            model = 1 + Group + Surf_idx + Group*Surf_idx + random(sub) + I;
            slm = SurfStatLinMod(Y, model);
            slm = SurfStatT(slm, surf_idx.*Group.IIa - surf_idx.*Group.IIb);
            p = tcdf(slm.t, slm.df);
            p = min([p, 1-p])
            p_value_set_2_10 = [ p_value_set_2_10 p];
            
            if(visualization_linearmodel)
                CorticalDepthVisualization(feature_idx, FCD_IIa_temp, FCD_IIb_temp, FCDIIa_feature, FCDIIb_feature, [-2.5 10], 1, 0.05, p_type_IIa_vs_IIb_2_10/2, p_type_IIa, p_type_IIb, FDR_threshold_final_cortical_depth_profiling);
                if(save_file_linearmodel)
                    export_fig([OUTPATH '/8_cortical_depth_analysis_FLAIR_sub_intensity_typeII_a_vs_b_CI' ], '-m3', '-png', '-painters'); close(gcf);
                end
            end
            
        end
        
    end
    
    for control_profling = 1
        
        % 1.1 Whole group patterns
        for all_control_analysis = 1
            
            feature_profile_cont = [];
            for i = 1 : length(feature_set_idx)                
                if(i == 1)
                    feature_profile_cont = [ feature_profile_cont ones(size(mean_z_lesion_cont, 3), 1)*NaN ];                
                end
                cont_data_temp = reshape(mean(mean_z_lesion_cont(feature_set_idx{i},:, :), 2), size(mean(mean_z_lesion_cont(feature_set_idx{i},:, :), 2), 1), size(mean(mean_z_lesion_cont(feature_set_idx{i},:, :), 2), 3));
                feature_profile_cont = [ feature_profile_cont mean(cont_data_temp, 1)' ];                                
                
            end
            
            label_feature_profile = {'Int', 'P-grad', 'T-grad', 'Sub Int', '', 'Thick', 'Curv', 'Depth', '', '', '', 'Int', 'P-grad', 'T-grad', 'Sub Int'};
            
            if(visualization)
                figure; set(gcf, 'Position', [587 500 1200 500]);
                hold on;
                linecolor = [0.8 0.8 0.8];
                hp1 = plot([0 52], [0 0]);
                hp2 = plot([0 52], [2 2], '--');
                hp3 = plot([0 52], [-2 -2], '--');
                set(hp1, 'Color', linecolor);
                set(hp2, 'Color', linecolor);
                set(hp3, 'Color', linecolor);
                set(gcf, 'Color', 'w');
                xtick = [0 3 6 9 12 16 19 22 27 30 33 36];
                h1 = notBoxPlot(feature_profile_cont, xtick, 0.5, 'line');
                flag_mu  = 'off';
                flag_sd  = 'off';
                flag_sem = 'off';
                barhead_linewidth = 0.3;
                for i = 2 : size(h1, 2)
                    set(h1(i).data, 'markerfacecolor',[1,1,1],'color',[0,0,0]);
                    set(h1(i).mu,   'visible' , flag_mu, 'markerfacecolor',[0,0,0],'color',[0,0,0]);
                    set(h1(i).sd,   'visible' , flag_sd, 'LineStyle', ':', 'Color',[0.3,0.3,0.3]);
                    set(h1(i).sem,  'visible' , flag_sem, 'LineStyle', '-', 'Color',[0.3,0.3,0.3]);
                    
                    x_1_temp = get(h1(i).mu, 'xData');
                    mu_1 = get(h1(i).mu, 'YData');
                    sd_1 = get(h1(i).sd, 'YData');
                    med_1 = median(get(h1(i).data, 'YData'));
                    
                    hold on; line([x_1_temp-barhead_linewidth x_1_temp+barhead_linewidth], [sd_1(1) sd_1(1)], 'LineWidth',2,'Color',[0 0 0]);
                    hold on; line([x_1_temp-barhead_linewidth x_1_temp+barhead_linewidth], [sd_1(2) sd_1(2)], 'LineWidth',2,'Color',[0 0 0]);
                    hold on; line([x_1_temp-barhead_linewidth x_1_temp+barhead_linewidth], [mu_1 mu_1], 'LineWidth',2,'Color',[0 0 0]);
                    hold on; scatter(x_1_temp, mu_1, 60, 'MarkerEdgeColor', [0 0 0], 'markerfacecolor', [0 0 0]);
                end
                
                xlim([2 max(xtick)+2]);
                set(gca, 'XTick', xtick);
                %% 0    3       6         9         12        16       19      22       27      30        33        36
                label_feature_profile_type = {'', 'Int', 'P-grad', 'T-grad', 'Sub Int', 'Thick', 'Curv', 'Depth', 'Int', 'P-grad', 'T-grad', 'Sub Int'};
                set(gca, 'XTickLabel', label_feature_profile_type);
                rotateticklabel(gca, 90);
                ylabel('z scores');
                ylim([-2.5 2.5]);
                set(gca, 'YTick', [-2.5 -2 0 2 2.5]);                
               
            end
            
        end        
        
        % 2 cortical depth analysis
        % 2-1. T1 intracortical intensity
        for depth_anlaysis_T1_cortical_intensity = 1
            
            excluded_case = { '306_1', '313_1', '322_1' };
            included_cases = find(~ismember(case_num_cont, excluded_case));
            
            feature_idx = 1;            
            feature_profile_cont= reshape(mean(mean_z_lesion_cont(feature_set_idx{feature_idx},:, :), 2), size(mean(mean_z_lesion_cont(feature_set_idx{feature_idx},:, :), 2), 1), size(mean(mean_z_lesion_cont(feature_set_idx{feature_idx},:, :), 2), 3));            
            
            if(visualization)
                x_lim = [0 3+(length(feature_set_idx{feature_idx})-1)*4];
                y_lim = [-4 4];
                
                visualization_graph_cortical_depth_profile_control(feature_set_idx{feature_idx}, feature_profile_cont', x_lim, y_lim);
                
            end            
         
        end
        
        % 2-2. T1 perpendicular gradient
        for depth_anlaysis_T1_perpendicular_gradient = 1
            
            excluded_case = { '306_1', '313_1', '322_1' };
            included_cases = find(~ismember(case_num_cont, excluded_case));
            
            feature_idx = 2;
            feature_profile_cont= reshape(mean(mean_z_lesion_cont(feature_set_idx{feature_idx},:, :), 2), size(mean(mean_z_lesion_cont(feature_set_idx{feature_idx},:, :), 2), 1), size(mean(mean_z_lesion_cont(feature_set_idx{feature_idx},:, :), 2), 3));
            
            if(visualization)
                x_lim = [0 3+(length(feature_set_idx{feature_idx})-1)*4];
                y_lim = [-4 4];
                
                visualization_graph_cortical_depth_profile_control(feature_set_idx{feature_idx}, feature_profile_cont', x_lim, y_lim);
                
            end
            
        end
        
        % 2-3. T1 tangential gradient
        for depth_anlaysis_T1_tangential_gradient = 1
            
            excluded_case = { '306_1', '313_1', '322_1' };
            included_cases = find(~ismember(case_num_cont, excluded_case));
            
            feature_idx = 3;
            feature_profile_cont= reshape(mean(mean_z_lesion_cont(feature_set_idx{feature_idx},:, :), 2), size(mean(mean_z_lesion_cont(feature_set_idx{feature_idx},:, :), 2), 1), size(mean(mean_z_lesion_cont(feature_set_idx{feature_idx},:, :), 2), 3));
            
            if(visualization)
                x_lim = [0 3+(length(feature_set_idx{feature_idx})-1)*4];
                y_lim = [-4 4];
                
                visualization_graph_cortical_depth_profile_control(feature_set_idx{feature_idx}, feature_profile_cont', x_lim, y_lim);
                
            end
            
        end
        
        % 2-4. T1 subcortical intensity
        for depth_anlaysis_T1_subcortical_intensity = 1
            
            excluded_case = { '306_1', '313_1', '322_1' };
            included_cases = find(~ismember(case_num_cont, excluded_case));
            
            feature_idx = 4;            
            feature_profile_cont= reshape(mean(mean_z_lesion_cont(feature_set_idx{feature_idx},:, :), 2), size(mean(mean_z_lesion_cont(feature_set_idx{feature_idx},:, :), 2), 1), size(mean(mean_z_lesion_cont(feature_set_idx{feature_idx},:, :), 2), 3));            
            
            if(visualization)
                x_lim = [0 3+(length(feature_set_idx{feature_idx})-1)*4];
                y_lim = [-4 4];
                
                visualization_graph_cortical_depth_profile_control(feature_set_idx{feature_idx}, feature_profile_cont', x_lim, y_lim);
                
            end     
            
        end
        
        % 2-5. FLAIR intracortical intensity
        for depth_anlaysis_FLAIR_cortical_intensity = 1
            
            excluded_case = { '306_1', '313_1', '322_1' };
            included_cases = find(~ismember(case_num_cont, excluded_case));
            
            feature_idx = 8;            
            feature_profile_cont= reshape(mean(mean_z_lesion_cont(feature_set_idx{feature_idx},:, :), 2), size(mean(mean_z_lesion_cont(feature_set_idx{feature_idx},:, :), 2), 1), size(mean(mean_z_lesion_cont(feature_set_idx{feature_idx},:, :), 2), 3));            
            
            if(visualization)
                x_lim = [0 3+(length(feature_set_idx{feature_idx})-1)*4];
                y_lim = [-4 4];
                
                visualization_graph_cortical_depth_profile_control(feature_set_idx{feature_idx}, feature_profile_cont', x_lim, y_lim);
                
            end     
            
        end
        
        % 2-6. FLAIR perpendicular gradient
        for depth_anlaysis_FLAIR_perpendicular_gradient = 1
            
            excluded_case = { '306_1', '313_1', '322_1' };
            included_cases = find(~ismember(case_num_cont, excluded_case));
            
            feature_idx = 9;            
            feature_profile_cont= reshape(mean(mean_z_lesion_cont(feature_set_idx{feature_idx},:, :), 2), size(mean(mean_z_lesion_cont(feature_set_idx{feature_idx},:, :), 2), 1), size(mean(mean_z_lesion_cont(feature_set_idx{feature_idx},:, :), 2), 3));            
            
            if(visualization)
                x_lim = [0 3+(length(feature_set_idx{feature_idx})-1)*4];
                y_lim = [-4 4];
                
                visualization_graph_cortical_depth_profile_control(feature_set_idx{feature_idx}, feature_profile_cont', x_lim, y_lim);
                
            end     
            
        end
        
        % 2-7. FLAIR tangential gradient
        for depth_anlaysis_FLAIR_tangential_gradient = 1
            
            excluded_case = { '306_1', '313_1', '322_1' };
            included_cases = find(~ismember(case_num_cont, excluded_case));
            
            feature_idx = 10;            
            feature_profile_cont= reshape(mean(mean_z_lesion_cont(feature_set_idx{feature_idx},:, :), 2), size(mean(mean_z_lesion_cont(feature_set_idx{feature_idx},:, :), 2), 1), size(mean(mean_z_lesion_cont(feature_set_idx{feature_idx},:, :), 2), 3));            
            
            if(visualization)
                x_lim = [0 3+(length(feature_set_idx{feature_idx})-1)*4];
                y_lim = [-4 4];
                
                visualization_graph_cortical_depth_profile_control(feature_set_idx{feature_idx}, feature_profile_cont', x_lim, y_lim);
                
            end     
            
        end
        
        % 2-8. FLAIR subcortical intensity
        for depth_anlaysis_FLAIR_subcortical_intensity = 1
            
            excluded_case = { '306_1', '313_1', '322_1' };
            included_cases = find(~ismember(case_num_cont, excluded_case));
            
            feature_idx = 11;            
            feature_profile_cont= reshape(mean(mean_z_lesion_cont(feature_set_idx{feature_idx},:, :), 2), size(mean(mean_z_lesion_cont(feature_set_idx{feature_idx},:, :), 2), 1), size(mean(mean_z_lesion_cont(feature_set_idx{feature_idx},:, :), 2), 3));            
            
            if(visualization)
                x_lim = [0 3+(length(feature_set_idx{feature_idx})-1)*4];
                y_lim = [-4 4];
                
                visualization_graph_cortical_depth_profile_control(feature_set_idx{feature_idx}, feature_profile_cont', x_lim, y_lim);
                
            end
            
        end
        
    end
    
end