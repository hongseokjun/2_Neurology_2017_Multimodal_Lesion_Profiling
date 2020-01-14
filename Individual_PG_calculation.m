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
end

for variable_dir_setup = 1
    process_dir = '/local_raid/seokjun/01_project/04_IntracorticalAnalysis/03_result/';
    outdir = '/local_raid/seokjun/01_project/04_IntracorticalAnalysis/03_result/validation_pg/';
   
    
    nummesh = 81920;
    intracortical_surface_number = 3;
    subcortical_surface_number = 3;
    sampling_space = 'native';
    pve_correction = 1;
    contrast_str = {'t1', 'flair'};
end

for process_pg = 1
    
    case_prefix = 'mcd';
    for i = 1 : size(case_num_pat, 1)
        case_num = case_num_pat{i};
        mkdir([ outdir '/' case_num ]);        
        for  j = 1 : 2
            if j == 1
                left_right = 'left';
            else
                left_right = 'right';
            end
            intracort_dir = '/local_raid/seokjun/01_project/04_IntracorticalAnalysis/03_result/';

            PerpendicularGradientCalculator_absolute_distance(process_dir, [ outdir '/' case_num ], case_prefix, case_num, left_right, nummesh, intracortical_surface_number, subcortical_surface_number, sampling_space, contrast_str);
            [r, s] = system(['DiffusionSmoother_absolute_distance.sh ' intracort_dir ' ' case_prefix ' ' case_num ' ' left_right ' ' num2str(intracortical_surface_number) ' ' num2str(subcortical_surface_number) ' 81920 ' sampling_space ' 5 quadratic 1 ' outdir ]);
            [r, s] = system(['SurfaceBasedResampler_absolute_distance.sh ' intracort_dir ' ' case_prefix ' ' case_num ' ' left_right ' ' num2str(intracortical_surface_number) ' ' num2str(subcortical_surface_number) ' 81920 ' sampling_space ' 5 quadratic ' outdir ]);
        end
    end
    
    case_prefix = 'TLE';
    for i = 19 : size(case_num_cont, 1)
        case_num = case_num_cont{i}
        mkdir([ outdir '/' case_num ]);
        for  j = 1 : 2
            if j == 1
                left_right = 'left';
            else
                left_right = 'right';
            end
            intracort_dir = '/local_raid/seokjun/01_project/04_IntracorticalAnalysis/03_result/';
            
            PerpendicularGradientCalculator_absolute_distance(process_dir, [ outdir '/' case_num ], case_prefix, case_num, left_right, nummesh, intracortical_surface_number, subcortical_surface_number, sampling_space, contrast_str);
            [r, s] = system(['DiffusionSmoother_absolute_distance.sh ' intracort_dir ' ' case_prefix ' ' case_num ' ' left_right ' ' num2str(intracortical_surface_number) ' ' num2str(subcortical_surface_number) ' 81920 ' sampling_space ' 5 quadratic 1 ' outdir ]);
            [r, s] = system(['SurfaceBasedResampler_absolute_distance.sh ' intracort_dir ' ' case_prefix ' ' case_num ' ' left_right ' ' num2str(intracortical_surface_number) ' ' num2str(subcortical_surface_number) ' 81920 ' sampling_space ' 5 quadratic ' outdir ]);
        end
    end
end