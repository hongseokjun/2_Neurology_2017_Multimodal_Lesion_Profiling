function PartialVolumeEffectRemover(output_dir, case_prefix, case_num, sampling_space, left_right, nummesh, intracortical_surface_number, contrast_str, csf_threshold_1, csf_threshold_2, ROI_kernel_size)

if(nummesh == 81920)
    numvert = 40962;
elseif(nummesh == 327680)
    numvert = 163842;
end

if(strcmp(sampling_space, 'native'))
    space='native';
    postfix_surf='_native';
    postfix_MRI='_nuc';
elseif(strcmp(sampling_space, 'tal'))
    space='final';
    postfix_surf='';
    postfix_MRI='_final';
end

for j = 1 : size(contrast_str, 2)
    postfix_contrast = [ '_' contrast_str{j} ];

     if(exist([ output_dir case_num '/measurement/' case_prefix '_' case_num '_intracortical_surface_' num2str(intracortical_surface_number) '_' left_right '_' num2str(nummesh) postfix_surf  postfix_contrast '_RI_corrected.txt' ], 'file') ~= 0)
          disp(['Already processed: ' case_prefix '_' case_num '_surfaces_' left_right '_' num2str(nummesh) postfix_surf  postfix_contrast '_RI_corrected.txt']);
          continue;
     end

    %% find WM peak ...
    filename=[ output_dir '/' case_num '/' space '/' case_prefix '_' case_num postfix_contrast '_WM_' space '.mnc.gz'];
    comp_flag = strfind(filename, '.mnc.gz');
    if(~isempty(comp_flag))
        [r, s] = system(['gunzip -vf ' filename]);
        filename = [ filename(1:comp_flag-1) ];
        vol = SurfStatReadVol1([filename '.mnc']);
    else
        vol = SurfStatReadVol1(filename);
    end

    [n,xout] = hist(vol.data(vol.data>=1 & vol.data<=99.5), 2000);
    n_sm = smooth(n);
    [peak_freq, peak_idx] = max(n_sm);
    WM_peak = xout(peak_idx)
    if(~isempty(comp_flag))
        [r, s] = system(['gzip -vf ' filename '.mnc']);
    end

    %% find GM peak ...
    filename=[ output_dir '/' case_num '/' space '/' case_prefix '_' case_num postfix_contrast '_GM_' space '.mnc.gz'];
    comp_flag = strfind(filename, '.mnc.gz');
    if(~isempty(comp_flag))
        [r, s] = system(['gunzip -vf ' filename]);
        filename = [ filename(1:comp_flag-1) ];
        vol = SurfStatReadVol1([filename '.mnc']);
    else
        vol = SurfStatReadVol1(filename);
    end

    [n,xout] = hist(vol.data(vol.data>=1 & vol.data<=99.5), 2000);
    n_sm = smooth(n);
    [peak_freq, peak_idx] = max(n_sm);
    GM_peak = xout(peak_idx)
    if(~isempty(comp_flag))
        [r, s] = system(['gzip -vf ' filename '.mnc']);
    end

    BG = (WM_peak + GM_peak)/2

    mri_name	   = [ output_dir '/' case_num '/' space '/' case_prefix '_' case_num postfix_contrast postfix_MRI '.mnc' ];
    pve_csf_name   = [ output_dir '/' case_num '/temp/' case_prefix '_' case_num postfix_contrast postfix_surf '_pve_csf.mnc' ];
    pve_gm_name    = [ output_dir '/' case_num '/temp/' case_prefix '_' case_num postfix_contrast postfix_surf '_pve_gm.mnc' ];
    csf_skel_name  = [ output_dir '/' case_num '/temp/' case_prefix '_' case_num postfix_contrast postfix_surf '_csf_skel.mnc' ];

    [r, s] = system(['gunzip -vf ' mri_name '.gz']);
    [r, s] = system(['gunzip -vf ' pve_csf_name '.gz']);
    [r, s] = system(['gunzip -vf ' pve_gm_name '.gz']);
    [r, s] = system(['gunzip -vf ' csf_skel_name '.gz']);

    MRI = SurfStatReadVol1(mri_name);
    pve_csf = SurfStatReadVol1(pve_csf_name);
    pve_gm  = SurfStatReadVol1(pve_gm_name);
    csf_skel = SurfStatReadVol1(csf_skel_name);

    [r, s] = system(['gzip -vf ' mri_name ]);
    [r, s] = system(['gzip -vf ' pve_csf_name ]);
    [r, s] = system(['gzip -vf ' pve_gm_name ]);
    [r, s] = system(['gzip -vf ' csf_skel_name ]);

    csf_skel_signal_vol = csf_skel.data * 0;
    csf_skel_signal_vol(csf_skel.data == 1 & pve_csf.data == 1 ) = MRI.data(csf_skel.data == 1 & pve_csf.data == 1);

    numcsf = zeros(intracortical_surface_number+2, numvert);
    cortical_surf_pve_csf = zeros(intracortical_surface_number+2, numvert);
    cortical_surf_pve_gm = zeros(intracortical_surface_number+2, numvert);

    for i = 1 : intracortical_surface_number + 2
        if(i == 1)
            basename  = [ case_prefix '_' case_num '_gray_surface_' left_right '_' num2str(nummesh) postfix_surf ];
        elseif(i == intracortical_surface_number + 2)
            basename  = [ case_prefix '_' case_num '_white_surface_' left_right '_' num2str(nummesh) postfix_surf ];
        else
            basename  = [ case_prefix '_' case_num '_intracortical_surface_' num2str(i-1) '_' left_right '_' num2str(nummesh) postfix_surf ];
        end

        cortical_surf_file{i} = [ output_dir '/' case_num '/surfaces/'  basename postfix_contrast '.obj' ];
        cortical_surf{i} = SurfStatReadSurf(cortical_surf_file{i});
        Edg=SurfStatEdg( cortical_surf{i} );
        
        [r, s] = system(['volume_object_evaluate -linear ' pve_csf_name ' ' cortical_surf_file{i} ' /tmp/tmp_' basename postfix_contrast '_pve_csf.txt' ]);
        [r, s] = system(['volume_object_evaluate -linear ' pve_gm_name  ' ' cortical_surf_file{i} ' /tmp/tmp_' basename postfix_contrast '_pve_gm.txt' ]);
        cortical_surf_pve_csf(i, :) = SurfStatReadData1(['/tmp/tmp_' basename postfix_contrast '_pve_csf.txt']);
        cortical_surf_pve_gm(i, :) = SurfStatReadData1(['/tmp/tmp_' basename postfix_contrast '_pve_gm.txt']);
        [r, s] = system(['rm -rf /tmp/tmp_' basename postfix_contrast '_pve_csf.txt']);
        [r, s] = system(['rm -rf /tmp/tmp_' basename postfix_contrast '_pve_gm.txt']);
    
        RI_values = SurfStatReadData1([ output_dir '/' case_num '/measurement/' basename postfix_contrast '_RI.txt' ]);

        if(i == 1)
            overcsf_index = find(cortical_surf_pve_csf(i, :) > csf_threshold_1);
            w_pure_gm     = 0.5;
            w_csf_average = 0.5;
        else
	    overcsf_index = find(cortical_surf_pve_csf(i, :) > csf_threshold_2);
            w_pure_gm     = 1;
            w_csf_average = 0;
        end

        [r, s] = system(['volume_object_evaluate -linear ' mri_name ' ' cortical_surf_file{i} ' /tmp/tmp_' basename postfix_contrast '_Intensity.txt' ]);
        image_signal = SurfStatReadData1(['/tmp/tmp_' basename postfix_contrast '_Intensity.txt']);
        [r, s] = system(['rm -rf /tmp/tmp_' basename postfix_contrast '_Intensity.txt']);
        image_signal_output = image_signal;
        size(overcsf_index, 2)
        for o = 1 : size(overcsf_index, 2)
            flag = 1;
            vox=(cortical_surf{i}.coord(:, overcsf_index(o))'-MRI.origin)./(MRI.vox(1:3))+1;
            voxel_coord = round(vox);
            ROI_kernel_size_temp = ROI_kernel_size;

            while(flag)
                x_start = (voxel_coord(1)-ROI_kernel_size_temp/2)*double(voxel_coord(1)-ROI_kernel_size_temp/2 > 1) + double(voxel_coord(1)-ROI_kernel_size_temp/2 <= 1);
                y_start = (voxel_coord(2)-ROI_kernel_size_temp/2)*double(voxel_coord(2)-ROI_kernel_size_temp/2 > 1) + double(voxel_coord(2)-ROI_kernel_size_temp/2 <= 1);
                z_start = (voxel_coord(3)-ROI_kernel_size_temp/2)*double(voxel_coord(3)-ROI_kernel_size_temp/2 > 1) + double(voxel_coord(3)-ROI_kernel_size_temp/2 <= 1);

                x_end = (x_start+ROI_kernel_size_temp)*double(x_start+ROI_kernel_size_temp <= MRI.dim(1)) + MRI.dim(1)*double(x_start+ROI_kernel_size_temp > MRI.dim(1));
                y_end = (y_start+ROI_kernel_size_temp)*double(y_start+ROI_kernel_size_temp <= MRI.dim(2)) + MRI.dim(2)*double(y_start+ROI_kernel_size_temp > MRI.dim(2));
                z_end = (z_start+ROI_kernel_size_temp)*double(z_start+ROI_kernel_size_temp <= MRI.dim(3)) + MRI.dim(3)*double(z_start+ROI_kernel_size_temp > MRI.dim(3));

                signal_kernel = csf_skel_signal_vol(int32(x_start:x_end), int32(y_start:y_end), int32(z_start:z_end));
                if(sum(signal_kernel(:) ~= 0) > 30)
                    csf_average = mean(signal_kernel(signal_kernel(:) ~= 0));
                    numcsf(i, overcsf_index(o)) = sum(signal_kernel(:) ~= 0);
                    flag = 0;
                else
                    ROI_kernel_size_temp = ROI_kernel_size_temp + 1;
                    flag = flag + 1;
                end
            end

            if(cortical_surf_pve_gm(i, overcsf_index(o)) == 0)
                pure_gm = 0;
            else
                if((1 - cortical_surf_pve_csf(i, overcsf_index(o)) - cortical_surf_pve_gm(i, overcsf_index(o))) < 0.3)
                    pure_gm = (image_signal(overcsf_index(o)) - cortical_surf_pve_csf(i, overcsf_index(o)) * csf_average)/cortical_surf_pve_gm(i, overcsf_index(o));
                else
                    pure_gm = (image_signal(overcsf_index(o)) - cortical_surf_pve_csf(i, overcsf_index(o)) * csf_average)/(1-cortical_surf_pve_csf(i, overcsf_index(o)));
                end
            end

            new_signal = w_pure_gm*pure_gm + w_csf_average*csf_average;
            if(cortical_surf_pve_csf(i, overcsf_index(o)) > 0.90)
                Ni=find_vertex_multineighbors2(cortical_surf{i}, overcsf_index(o), 3, Edg);
                new_signal = mean([image_signal_output(Ni) new_signal]);
            end

            if(new_signal > 100)
                new_signal = 100;
            end
            if(new_signal < 0)
                new_signal = 0;
            end

            image_signal_output(overcsf_index(o)) = new_signal;
            RI_values(overcsf_index(o)) = 100.0 * (new_signal - GM_peak) / abs(BG - GM_peak);
        end

        SurfStatWriteData( [ output_dir '/' case_num '/measurement/' basename postfix_contrast '_RI_corrected.txt' ], RI_values, 'a');
        SurfStatWriteData( [ output_dir '/' case_num '/measurement/' basename postfix_contrast '_pve_csf.txt' ], cortical_surf_pve_csf(i, :), 'a');
    end
end
