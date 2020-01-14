clear all;
close all;

PREFIX    = 'mcd';
datadir   = '/local_raid/seokjun/01_project/04_IntracorticalAnalysis/03_result/inter_rater_reliability/data/final_data/';
Cases_pat = '/local_raid/seokjun/01_project/04_IntracorticalAnalysis/01_analysis/misc/caselist_IRR2.txt'; % caselist_IRR1.txt or caselist_IRR2.txt

fid = fopen(Cases_pat);
demo = textscan(fid, '%s', 'Delimiter', ',', 'CollectOutput', 1);
case_num_pat = demo{1};

vol_DS  = zeros(181, 217, 181, length(case_num_pat));
vol_NAB = zeros(181, 217, 181, length(case_num_pat));
for i = 1 : length(case_num_pat)
    
    i
    
    case_curr = case_num_pat{i};
    vol = SurfStatReadVol1([ datadir '/' PREFIX '_' case_curr '_DS_new.mnc' ]);   vol_DS(:, :, :, i) = vol.data;
    vol = SurfStatReadVol1([ datadir '/' PREFIX '_' case_curr '_NAB_new.mnc' ]);  vol_NAB(:, :, :, i) = vol.data;
    
end

diceidx = zeros(length(case_num_pat), 1);
for i = 1 : length(case_num_pat)
    
    i
    
    vol_DS_temp  = vol_DS(:,:,:, i);
    vol_NAB_temp = vol_NAB(:,:,:, i);
    
    vol_DS_temp  = vol_DS_temp(:);
    vol_NAB_temp = vol_NAB_temp(:);
    
    M1 = sum(vol_DS_temp);
    M2 = sum(vol_NAB_temp);
    
    M1_M2 = sum(vol_DS_temp.*vol_NAB_temp);
    
    diceidx(i) = 2*M1_M2/(M1+M2);
    
end

disp([ 'N=' num2str(length(diceidx)) ', mean±std(range)=' num2str(mean(diceidx)) '±' num2str(std(diceidx)) '(' num2str(min(diceidx)) '-' num2str(max(diceidx)) ')' ]);