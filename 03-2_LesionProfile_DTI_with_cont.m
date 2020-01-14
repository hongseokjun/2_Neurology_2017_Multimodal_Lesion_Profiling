clear all
close all
addpath('/data/noel/noel2/local/matlab/surfstat_chicago');

%% The difference with the original script is that here, we have control data that were averaged across lesion locations.
%% Therefore we are able to run non-parametric tests. But you should note that this way to generate control data can artificially
%% amplify the statistical significances as it actually reduces the standard deviation of controls to less than 1.

%% Read pre-calculated zscore of T1 and FLAIR in the lesion
for read_matfiles = 1
    
    load('zscore_database_DTI.mat');
    load('zscore_database_DTI_cont.mat');
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
    Group_pat                 = 'FCD'
    Prefix_pat                = 'mcd'
    Cases_cont                = '/local_raid/seokjun/01_project/04_IntracorticalAnalysis/01_analysis/Demographic_data_control_DTI.txt'
    Cases_pat                 = '/local_raid/seokjun/01_project/04_IntracorticalAnalysis/01_analysis/Demographic_data_FCD_DTI.txt'
    
    Left_Right                = 'both'
    NumIntSurf                = 3
    NumSubSurf                = 4
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
    
    OUTPATH = '/local_raid/seokjun/01_project/04_IntracorticalAnalysis/03_result/lesion_profile/';
    case_num_pat_temp = case_num_pat;
    case_cont_real = case_num_cont';
    mean_z_lesion_org = mean_z_lesion;
    
    FDR_threshold_final_lesion_profiling = FDR([p_value_lesion_profile_T1_FLAIR p_value_lesion_profile_DTI], 0.05); %% 0.0114
    FDR_threshold_final_cortical_depth_profiling = FDR([p_value_lesion_cortical_depth_profile p_value_lesion_cortical_depth_profile_DTI], 0.05); %% 0.02
    
end

%% Select features to be included in the analyses
for feature_selection = 1
    
    %% mean_z_lesion
    %% 1:5   = FA of DTI Intra         -> 3 5       p1
    %% 6:12  = FA of DTI Subcortical   -> 7 9 11    p2
    %% 13:17 = MD of DTI Intra         -> 12        p3
    %% 18:24 = MD of DTI Subcortical   -> 14 16 18  p4
    feature_set_idx = { 3, 5, 7, 9, 11, 15, 17, 19, 21, 23 };
    include = ones(size(case_num_pat_temp, 1), 1);
    
end

%% Start profile analyses
for profile_analyses = 1
 
    for patient_profiling = 1
        
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
            
            p_value_set_1_1 = p_type_II;
            
            feature_profile = [
                ones(1, length(case_num_pat_temp))*NaN
                mean(mean_z_lesion(3, find(include ~= 0)), 1)
                mean(mean_z_lesion(5, find(include ~= 0)), 1)
                mean(mean_z_lesion(7, find(include ~= 0)), 1)
                mean(mean_z_lesion(9, find(include ~= 0)), 1)
                mean(mean_z_lesion(11, find(include ~= 0)), 1)
                mean(mean_z_lesion(15, find(include ~= 0)), 1)
                mean(mean_z_lesion(17, find(include ~= 0)), 1)
                mean(mean_z_lesion(19, find(include ~= 0)), 1)
                mean(mean_z_lesion(21, find(include ~= 0)), 1)
                mean(mean_z_lesion(23, find(include ~= 0)), 1)
                ]';
            feature_profile(feature_profile > 10) = 3;
            
            if(visualization)
                label_feature_profile = {'mid', 'WM', '2mm', '4mm', '6mm' 'mid', 'WM', '2mm', '4mm', '6mm'};
                
                figure; set(gcf, 'Position', [587 600 1200 500]);
                hold on;
                linecolor = [0.8 0.8 0.8];
                hp1 = plot([0 52], [0 0]);
                hp2 = plot([0 52], [2 2], '--');
                hp3 = plot([0 52], [-2 -2], '--');
                set(hp1, 'Color', linecolor);
                set(hp2, 'Color', linecolor);
                set(hp3, 'Color', linecolor);
                set(gcf, 'Color', 'w');
                xtick = [0 3 6 9 12 15 18 21 24 27 30 ];
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
                % 0     3      6      9     12    15     18    21     24     27     30
                label_feature_profile_type = {'', 'mid', 'WM', '2mm', '4mm', '6mm' 'mid', 'WM', '2mm', '4mm', '6mm'};
                
                set(gca, 'XTickLabel', label_feature_profile_type);
                ylabel('z scores');
                ylim([-5 10]);
                set(gca, 'YTick', [-5 -2 0 2 5]);
                rotateticklabel(gca, 90);
                
                hold on;
                x_range = xtick(2:end);
                FDR_threshold = FDR_threshold_final_lesion_profiling;
                for i = 1 : size(p_value_set_1_1, 2)
                    if(p_value_set_1_1(i) <= FDR_threshold)
                        scatter(x_range(i), 9, 50, 'rv');
                    elseif(p_value_set_1_1(i) <= 0.05 && p_value_set_1_1(i) > FDR_threshold)
                        scatter(x_range(i), 9, 50, 'r*');
                    end
                end
                
                if(save_file)
                    export_fig([OUTPATH '/9_cortical_depth_analysis_DTI_FA_MD' ], '-m3', '-png', '-painters'); close(gcf);
                end
            end
            
        end
        
        % 1.2 Group separation along the histological types: FCD IIa, IIb
        for separation_analysis_by_histo_types = 1
            
            FCD_IIa = strncmp(histo_type, 'FCDIIa', 6)';
            FCD_IIb = strncmp(histo_type, 'FCDIIb', 6)';
            
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
            p_value_set_1_2 = [p_type_IIa p_type_IIb];
            
            p_type_IIa_vs_IIb_1_2 = [];
            include = ones(size(case_num_pat_temp, 1), 1);
            for i = 1 : length(feature_set_idx)
                [p_temp,h,stats] = ranksum(mean(mean_z_lesion(feature_set_idx{i}, FCD_IIa), 1), mean(mean_z_lesion(feature_set_idx{i}, FCD_IIb), 1));
                p_type_IIa_vs_IIb_1_2 = [ p_type_IIa_vs_IIb_1_2 p_temp ];
            end
            
            feature_profile_typeIIa = [
                ones(1, size(find(FCD_IIa), 2))*NaN
                mean(mean_z_lesion(3, FCD_IIa), 1)
                mean(mean_z_lesion(5, FCD_IIa), 1)
                mean(mean_z_lesion(7, FCD_IIa), 1)
                mean(mean_z_lesion(9, FCD_IIa), 1)
                mean(mean_z_lesion(11, FCD_IIa), 1)
                mean(mean_z_lesion(15, FCD_IIa), 1)
                mean(mean_z_lesion(17, FCD_IIa), 1)
                mean(mean_z_lesion(19, FCD_IIa), 1)
                mean(mean_z_lesion(21, FCD_IIa), 1)
                mean(mean_z_lesion(23, FCD_IIa), 1)
                ]';
            feature_profile_typeIIa(feature_profile_typeIIa > 5) = 3;
            
            feature_profile_typeIIb = [
                ones(1, size(find(FCD_IIb), 2))*NaN
                mean(mean_z_lesion(3, FCD_IIb), 1)
                mean(mean_z_lesion(5, FCD_IIb), 1)
                mean(mean_z_lesion(7, FCD_IIb), 1)
                mean(mean_z_lesion(9, FCD_IIb), 1)
                mean(mean_z_lesion(11, FCD_IIb), 1)
                mean(mean_z_lesion(15, FCD_IIb), 1)
                mean(mean_z_lesion(17, FCD_IIb), 1)
                mean(mean_z_lesion(19, FCD_IIb), 1)
                mean(mean_z_lesion(21, FCD_IIb), 1)
                mean(mean_z_lesion(23, FCD_IIb), 1)
                ]';
            
            feature_profile_typeIIb(feature_profile_typeIIb > 5) = 3;
            
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
                xtick = [0 3 6 9 12 15 18 21 24 27 30 ];
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
                
                % 0     3      6      9     12    15     18    21     24     27     30
                label_feature_profile_type = {'', 'mid', 'WM', '2mm', '4mm', '6mm' 'mid', 'WM', '2mm', '4mm', '6mm'};
                ylabel('z scores');
                ylim([-5 11]);
                set(gca, 'XTickLabel', label_feature_profile_type);
                set(gca, 'YTick', [ -5 -2 0 2 5 ]);
                rotateticklabel(gca, 90);
                
                hold on;
                x_range = xtick(2:end);
                FDR_threshold = FDR_threshold_final_cortical_depth_profiling;
                for i = 1 : size(p_type_IIa, 2)
                    if(p_type_IIa(i) <= FDR_threshold)
                        scatter(x_range(i), -4, 50, 'k^');
                    elseif(p_type_IIa(i) < 0.05 && p_type_IIa(i) > FDR_threshold)
                        scatter(x_range(i), -4, 50, 'k*');
                    end
                end
                
                for i = 1 : size(p_type_IIb, 2)
                    if(p_type_IIb(i) <= FDR_threshold)
                        scatter(x_range(i)+1, -4, 50, 'k^');
                    elseif(p_type_IIb(i) < 0.05 && p_type_IIb(i) > FDR_threshold)
                        scatter(x_range(i)+1, -4, 50, 'k*');
                    end
                end
                
                if(save_file)
                    export_fig([OUTPATH '/10_cortical_depth_analysis_DTI_FA_MD_typeIIa_IIb' ], '-m3', '-png', '-painters'); close(gcf);
                end
            end
            
            p_value_set_1_2 = [p_value_set_1_2 p_type_IIa_vs_IIb_1_2];
            
        end
        
        % 1.3 Group separation along the histological types: FCD IIa, IIb
        for GM_WM_wise_separation_analysis_by_histo_types = 1
            
            feature_set_idx = { [3 5], [7, 9], [15 17], [19, 21] };
            FCD_IIa = strncmp(histo_type, 'FCDIIa', 6)';
            FCD_IIb = strncmp(histo_type, 'FCDIIb', 6)';
            
            p_type_IIa = [];
            p_type_IIb = [];
            for i = 1 : length(feature_set_idx)
                pat_data  = mean(mean_z_lesion(feature_set_idx{i}, FCD_IIa), 1);
                cont_data = mean(reshape(mean(mean_z_lesion_cont(feature_set_idx{i}, :, :), 2), size(mean(mean_z_lesion_cont(feature_set_idx{i}, :, :), 2), 1), size(mean(mean_z_lesion_cont(feature_set_idx{i}, :, :), 2), 3)), 1);
                [ p h stats ] = ranksum(pat_data, cont_data);
                
                p_type_IIa = [ p_type_IIa p ];
                
                
                pat_data  = mean(mean_z_lesion(feature_set_idx{i}, FCD_IIb), 1);
                cont_data = mean(reshape(mean(mean_z_lesion_cont(feature_set_idx{i}, :, :), 2), size(mean(mean_z_lesion_cont(feature_set_idx{i}, :, :), 2), 1), size(mean(mean_z_lesion_cont(feature_set_idx{i}, :, :), 2), 3)), 1);
                [ p h stats ] = ranksum(pat_data, cont_data);
                
                p_type_IIb = [ p_type_IIb p ];
            end
            
            p_value_set_2_1 = [p_type_IIa p_type_IIb];
            
            include = ones(size(case_num_pat_temp, 1), 1);
            p_type_IIa_vs_IIb_2_1  = [];
            
            for i = 1 : length(feature_set_idx)
                pat_data_1  = mean(mean_z_lesion(feature_set_idx{i}, FCD_IIa), 1);
                pat_data_2  = mean(mean_z_lesion(feature_set_idx{i}, FCD_IIb), 1);
                
                [p_temp,h,stats] = ranksum(pat_data_1, pat_data_2);
                p_type_IIa_vs_IIb_2_1 = [ p_type_IIa_vs_IIb_2_1 p_temp ];
            end
            
            feature_profile_typeIIa = [
                ones(1, size(find(FCD_IIa), 2))*NaN
                mean(mean_z_lesion(feature_set_idx{1}, FCD_IIa), 1)
                mean(mean_z_lesion(feature_set_idx{2}, FCD_IIa), 1)
                mean(mean_z_lesion(feature_set_idx{3}, FCD_IIa), 1)
                mean(mean_z_lesion(feature_set_idx{4}, FCD_IIa), 1)
                ]';
            feature_profile_typeIIa(feature_profile_typeIIa > 5) = 3;
            
            feature_profile_typeIIb = [
                ones(1, size(find(FCD_IIb), 2))*NaN
                mean(mean_z_lesion(feature_set_idx{1}, FCD_IIb), 1)
                mean(mean_z_lesion(feature_set_idx{2}, FCD_IIb), 1)
                mean(mean_z_lesion(feature_set_idx{3}, FCD_IIb), 1)
                mean(mean_z_lesion(feature_set_idx{4}, FCD_IIb), 1)
                ]';
            
            feature_profile_typeIIb(feature_profile_typeIIb > 5) = 3;
            
            if(visualization)
                f = figure; set(gcf, 'Position', [587 500 700 500]);
                hold on;
                linecolor = [0.8 0.8 0.8];
                hp1 = plot([0 52], [0 0]);
                hp2 = plot([0 52], [2 2], '--');
                hp3 = plot([0 52], [-2 -2], '--');
                set(hp1, 'Color', linecolor);
                set(hp2, 'Color', linecolor);
                set(hp3, 'Color', linecolor);
                set(gcf, 'Color', 'w');
                xtick = [0 3 6 9 12 ];
                h1 = notBoxPlot(feature_profile_typeIIa, xtick, 0.5, 'line');
                h2 = notBoxPlot(feature_profile_typeIIb, xtick+1, 0.5, 'line');
                flag_mu  = 'off';
                flag_sd  = 'off';
                flag_sem = 'off';
                barhead_linewidth = 0.2;
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
                
                % 0     3    6      9     12
                label_feature_profile_type = {'', 'Int', 'WM', 'Int', 'WM' };
                ylabel('z scores');
                ylim([-7 11]);
                set(gca, 'XTickLabel', label_feature_profile_type);
                set(gca, 'YTick', [ -5 -2 0 2 5 ]);
                % rotateticklabel(gca, 90);
                
                hold on;
                x_range = xtick(2:end);
                FDR_threshold = FDR_threshold_final_lesion_profiling;
                for i = 1 : size(p_type_IIa, 2)
                    if(p_type_IIa(i) <= FDR_threshold)
                        scatter(x_range(i), -4, 50, 'k^');
                    elseif(p_type_IIa(i) < 0.05 && p_type_IIa(i) > FDR_threshold)
                        scatter(x_range(i), -4, 50, 'k*');
                    end
                end
                
                for i = 1 : size(p_type_IIb, 2)
                    if(p_type_IIb(i) <= FDR_threshold)
                        scatter(x_range(i)+1, -4, 50, 'k^');
                    elseif(p_type_IIb(i) < 0.05 && p_type_IIb(i) > FDR_threshold)
                        scatter(x_range(i)+1, -4, 50, 'k*');
                    end
                end
                
                if(save_file)
                    export_fig([OUTPATH '/10_FCD_lesion_feature_profile_DTI_FA_MD_typeIIa_IIb_new' ], '-m4', '-png' ); close(gcf);
                end
            end
            
            p_value_set_2_1 = [p_value_set_2_1 p_type_IIa_vs_IIb_2_1];
            
        end
        
        % 2 FCD IIa vs. IIb - cortical depth analysis
        % 2-1. FA
        for depth_anlaysis_FA = 1
            
            excluded_case = { };
            included_cases = ~ismember(case_num_pat_temp, excluded_case);
            FCD_IIa_temp = FCD_IIa' & included_cases;
            FCD_IIb_temp = FCD_IIb' & included_cases;
            
            feature_idx = [ 3 5 7 9 ] ;
            feature_profile_temp = mean_z_lesion(feature_idx, :);
            feature_profile_cont_temp = reshape(mean(mean_z_lesion_cont(feature_idx, :, :), 2), size(mean(mean_z_lesion_cont(feature_idx, :, :), 2), 1), size(mean(mean_z_lesion_cont(feature_idx, :, :), 2), 3));
            FCDIIa_feature = feature_profile_temp(:, FCD_IIa_temp)';
            FCDIIb_feature = feature_profile_temp(:, FCD_IIb_temp)';
            
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
                y_lim = [-3 4];
                
                visualization_graph_cortical_depth_profile(feature_idx, p_value_set_2_2, FCDIIa_feature, FCDIIb_feature, x_lim, y_lim, FDR_threshold_final_cortical_depth_profiling);
                
                if(save_file)
                    export_fig([OUTPATH '/12_cortical_depth_analysis_DTI_FA_typeII_a_vs_b_new2' ], '-m3', '-png', '-painters'); close(gcf);
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
                CorticalDepthVisualization(feature_idx, FCD_IIa_temp, FCD_IIb_temp, FCDIIa_feature, FCDIIb_feature, [-5 7], 1, 0.05, p_type_IIa_vs_IIb_2_1/2, p_type_IIa, p_type_IIb, FDR_threshold_final_cortical_depth_profiling);
                if(save_file_linearmodel)
                    export_fig([OUTPATH '/12_cortical_depth_analysis_FA_typeII_a_vs_b_CI' ], '-m3', '-png', '-painters'); close(gcf);
                end
            end
            
        end
        
        % 2-2. MD
        for depth_anlaysis_MD = 1
            
            excluded_case = { };
            included_cases = ~ismember(case_num_pat_temp, excluded_case);
            FCD_IIa_temp = FCD_IIa' & included_cases;
            FCD_IIb_temp = FCD_IIb' & included_cases;
            
            feature_idx = [ 15 17 19 21 ] ;
            feature_profile_temp = mean_z_lesion(feature_idx, :);
            feature_profile_temp(feature_profile_temp>5) = 5;
            feature_profile_cont_temp = reshape(mean(mean_z_lesion_cont(feature_idx, :, :), 2), size(mean(mean_z_lesion_cont(feature_idx, :, :), 2), 1), size(mean(mean_z_lesion_cont(feature_idx, :, :), 2), 3));
            FCDIIa_feature = feature_profile_temp(:, FCD_IIa_temp)';
            FCDIIb_feature = feature_profile_temp(:, FCD_IIb_temp)';
            
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
                y_lim = [-2 4];
                
                visualization_graph_cortical_depth_profile(feature_idx, p_value_set_2_3, FCDIIa_feature, FCDIIb_feature, x_lim, y_lim, FDR_threshold_final_cortical_depth_profiling);
                
                if(save_file)
                    export_fig([OUTPATH '/13_cortical_depth_analysis_DTI_MD_typeII_a_vs_b_new2' ], '-m3', '-png', '-painters'); close(gcf);
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
                CorticalDepthVisualization(feature_idx, FCD_IIa_temp, FCD_IIb_temp, FCDIIa_feature, FCDIIb_feature, [-5 7], 1, 0.05, p_type_IIa_vs_IIb_2_1/2, p_type_IIa, p_type_IIb, FDR_threshold_final_cortical_depth_profiling);
                if(save_file_linearmodel)
                    export_fig([OUTPATH '/13_cortical_depth_analysis_MD_typeII_a_vs_b_CI' ], '-m3', '-png', '-painters'); close(gcf);
                end
            end
            
        end
        
    end
    
    for control_profiling = 1
        
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
            
            label_feature_profile = {'mid', 'WM', '2mm', '4mm', '6mm' 'mid', 'WM', '2mm', '4mm', '6mm'};
            
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
                xtick = [0 3 6 9 12 15 18 21 24 27 30 ];
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
                % 0     3      6      9     12    15     18    21     24     27     30
                label_feature_profile_type = {'', 'mid', 'WM', '2mm', '4mm', '6mm' 'mid', 'WM', '2mm', '4mm', '6mm'};
                set(gca, 'XTickLabel', label_feature_profile_type);
                rotateticklabel(gca, 90);
                ylabel('z scores');
                ylim([-2.5 2.5]);
                set(gca, 'YTick', [-2.5 -2 0 2 2.5]);                
               
            end            
            
        end
        
               
        % 2 cortical depth analysis
        % 2-1. FA
        for depth_anlaysis_FA = 1
            
            excluded_case = { '306_1', '313_1', '322_1' };
            included_cases = find(~ismember(case_num_cont, excluded_case));
            
            feature_idx = [ 3 5 7 9 ] ;            
            feature_profile_cont= reshape(mean(mean_z_lesion_cont(feature_idx,:, :), 2), size(mean(mean_z_lesion_cont(feature_idx,:, :), 2), 1), size(mean(mean_z_lesion_cont(feature_idx,:, :), 2), 3));            
            
            if(visualization)
                x_lim = [0 3+(length(feature_idx)-1)*4];
                y_lim = [-4 4];
                
                visualization_graph_cortical_depth_profile_control(feature_idx, feature_profile_cont', x_lim, y_lim);
                
            end
            
        end
        
        % 2-2. MD
        for depth_anlaysis_MD = 1
            
            excluded_case = { '306_1', '313_1', '322_1' };
            included_cases = find(~ismember(case_num_cont, excluded_case));
            
            feature_idx = [ 15 17 19 21 ] ;            
            feature_profile_cont= reshape(mean(mean_z_lesion_cont(feature_idx,:, :), 2), size(mean(mean_z_lesion_cont(feature_idx,:, :), 2), 1), size(mean(mean_z_lesion_cont(feature_idx,:, :), 2), 3));            
            
            if(visualization)
                x_lim = [0 3+(length(feature_idx)-1)*4];
                y_lim = [-4 4];
                
                visualization_graph_cortical_depth_profile_control(feature_idx, feature_profile_cont', x_lim, y_lim);
                
            end

        end        
    end
    
end