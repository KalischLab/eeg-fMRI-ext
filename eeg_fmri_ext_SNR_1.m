%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% eeg_fmri_ext_SNR_1.m
%
% written by B.Meyer, Neuroimaging Center Mainz, 2026
%
% eeg_fmri_ext_SNR_1 BOLD realigns each echo's time series to remove motion
% artifacts before SNR calculation. T1 images are also segmented to gm-mask
% SNR files.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all

%Data is copied to this temporary directory as data cannot be processed on
%shared drive (nicshare)

realignment_directory_temp = '/home/mey1b/EEG_fmri_ext/data/temp/';
realignment_directory_out = '/home/mey1b/EEG_fmri_ext/data/EEG_fMRI_ext_realigned_for_tSNR/';

%restscans (1,6: without EEG, 2,3,4,5: with EEG)
restscans = 1:6;
%subjects
subjects = 10:34;

if ~exist(realignment_directory_out,'dir')
    mkdir(realignment_directory_out)
end

for sub = subjects
    for rest = restscans

        switch rest
            case 1
                rest_str = 'rest1';
                day = 1;
            case 2
                rest_str = 'rest2';
                day = 2;
            case 3
                rest_str = 'rest3';
                day = 2;
            case 4
                rest_str = 'rest4';
                day = 2;
            case 5
                rest_str = 'rest5';
                day = 2;
            case 6
                rest_str = 'rest6';
                day = 3;
        end

        for echo = 1:4
            if ~exist(realignment_directory_temp,'dir')
                mkdir(realignment_directory_temp)
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Copy and unzip the functional images
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            path_func               = ['/mnt/nicshare/EEG_fmri_ext/data/bids/Nifti/sub-S' num2str(sub) '/ses-Day' num2str(day) '/func/'];
            filename_func           = ['sub-S' num2str(sub) '_ses-Day' num2str(day) '_task-' rest_str '_acq-mb3me4_echo-' num2str(echo) '_bold.nii.gz'];
            filename_unzipped_func  = ['sub-S' num2str(sub) '_ses-Day' num2str(day) '_task-' rest_str '_acq-mb3me4_echo-' num2str(echo) '_bold.nii'];
            origin_func             = [path_func filename_func];

            if ~exist(origin_func,'file')
                continue
            end

            destination_func = [realignment_directory_temp filename_func];
            copyfile(origin_func,destination_func);
            gunzip(destination_func);
            destination_func = [realignment_directory_temp filename_unzipped_func];

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Copy and unzip the anatomical image
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            path_anat               = ['/mnt/nicshare/EEG_fmri_ext/data/bids/Nifti/sub-S' num2str(sub) '/ses-Day1/anat/'];
            filename_anat           = ['sub-S' num2str(sub) '_ses-Day1_T1w.nii.gz'];
            filename_unzipped_anat  = ['sub-S' num2str(sub) '_ses-Day1_T1w.nii'];
            origin_anat             = [path_anat filename_anat];

            if ~exist(origin_anat,'file')
                continue
            end

            destination_anat = [realignment_directory_out filename_anat];

            copyfile(origin_anat,destination_anat);
            gunzip(destination_anat);

            destination_anat = [realignment_directory_out filename_unzipped_anat];

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Run SPM batch script to perform realignment of all echos. When
            %first echo is processed a segmentation is performed to obtain
            %GM and WM masks.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


            if echo == 1
                matlabbatch = [];
                matlabbatch{1}.spm.util.split.vol = {[destination_func ',1']};
                matlabbatch{1}.spm.util.split.outdir = {realignment_directory_temp};
                matlabbatch{2}.spm.spatial.realign.estwrite.data{1}(1) = cfg_dep('4D to 3D File Conversion: Series of 3D Volumes', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','splitfiles'));
                matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
                matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.sep = 4;
                matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
                matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
                matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.interp = 2;
                matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
                matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.weight = '';
                matlabbatch{2}.spm.spatial.realign.estwrite.roptions.which = [0 1];
                matlabbatch{2}.spm.spatial.realign.estwrite.roptions.interp = 4;
                matlabbatch{2}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
                matlabbatch{2}.spm.spatial.realign.estwrite.roptions.mask = 1;
                matlabbatch{2}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
                matlabbatch{3}.spm.util.cat.vols(1) = cfg_dep('Realign: Estimate & Reslice: Realigned Images (Sess 1)', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{1}, '.','cfiles'));
                matlabbatch{3}.spm.util.cat.name = [realignment_directory_out '/sub-S' num2str(sub) '_ses-Day' num2str(day) '_task-' rest_str '_acq-mb3me4_echo-' num2str(echo) '_bold_realigned.nii'];
                matlabbatch{3}.spm.util.cat.dtype = 4;
                matlabbatch{3}.spm.util.cat.RT = NaN;
                matlabbatch{4}.spm.spatial.coreg.estimate.ref(1) = cfg_dep('Realign: Estimate & Reslice: Mean Image', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','rmean'));
                matlabbatch{4}.spm.spatial.coreg.estimate.source = {[destination_anat ',1']};
                matlabbatch{4}.spm.spatial.coreg.estimate.other = {''};
                matlabbatch{4}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
                matlabbatch{4}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
                matlabbatch{4}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
                matlabbatch{4}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];

                matlabbatch{5}.spm.spatial.preproc.channel.vols(1) = cfg_dep('Coregister: Estimate: Coregistered Images', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','cfiles'));
                matlabbatch{5}.spm.spatial.preproc.channel.biasreg = 0.001;
                matlabbatch{5}.spm.spatial.preproc.channel.biasfwhm = 60;
                matlabbatch{5}.spm.spatial.preproc.channel.write = [0 0];
           
                matlabbatch{5}.spm.spatial.preproc.tissue(1).tpm = {'/usr/local/MATLAB/Toolboxes/SPM12/spm/tpm/TPM.nii,1'};
                matlabbatch{5}.spm.spatial.preproc.tissue(1).ngaus = 1;
                matlabbatch{5}.spm.spatial.preproc.tissue(1).native = [1 0];
                matlabbatch{5}.spm.spatial.preproc.tissue(1).warped = [0 0];
                matlabbatch{5}.spm.spatial.preproc.tissue(2).tpm = {'/usr/local/MATLAB/Toolboxes/SPM12/spm/tpm/TPM.nii,2'};
                matlabbatch{5}.spm.spatial.preproc.tissue(2).ngaus = 1;
                matlabbatch{5}.spm.spatial.preproc.tissue(2).native = [1 0];
                matlabbatch{5}.spm.spatial.preproc.tissue(2).warped = [0 0];
                matlabbatch{5}.spm.spatial.preproc.tissue(3).tpm = {'/usr/local/MATLAB/Toolboxes/SPM12/spm/tpm/TPM.nii,3'};
                matlabbatch{5}.spm.spatial.preproc.tissue(3).ngaus = 2;
                matlabbatch{5}.spm.spatial.preproc.tissue(3).native = [1 0];
                matlabbatch{5}.spm.spatial.preproc.tissue(3).warped = [0 0];
                matlabbatch{5}.spm.spatial.preproc.tissue(4).tpm = {'/usr/local/MATLAB/Toolboxes/SPM12/spm/tpm/TPM.nii,4'};
                matlabbatch{5}.spm.spatial.preproc.tissue(4).ngaus = 3;
                matlabbatch{5}.spm.spatial.preproc.tissue(4).native = [1 0];
                matlabbatch{5}.spm.spatial.preproc.tissue(4).warped = [0 0];
                matlabbatch{5}.spm.spatial.preproc.tissue(5).tpm = {'/usr/local/MATLAB/Toolboxes/SPM12/spm/tpm/TPM.nii,5'};
                matlabbatch{5}.spm.spatial.preproc.tissue(5).ngaus = 4;
                matlabbatch{5}.spm.spatial.preproc.tissue(5).native = [1 0];
                matlabbatch{5}.spm.spatial.preproc.tissue(5).warped = [0 0];
                matlabbatch{5}.spm.spatial.preproc.tissue(6).tpm = {'/usr/local/MATLAB/Toolboxes/SPM12/spm/tpm/TPM.nii,6'};
                matlabbatch{5}.spm.spatial.preproc.tissue(6).ngaus = 2;
                matlabbatch{5}.spm.spatial.preproc.tissue(6).native = [0 0];
                matlabbatch{5}.spm.spatial.preproc.tissue(6).warped = [0 0];
                matlabbatch{5}.spm.spatial.preproc.warp.mrf = 1;
                matlabbatch{5}.spm.spatial.preproc.warp.cleanup = 1;
                matlabbatch{5}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
                matlabbatch{5}.spm.spatial.preproc.warp.affreg = 'mni';
                matlabbatch{5}.spm.spatial.preproc.warp.fwhm = 0;
                matlabbatch{5}.spm.spatial.preproc.warp.samp = 3;
                matlabbatch{5}.spm.spatial.preproc.warp.write = [0 0];
                matlabbatch{5}.spm.spatial.preproc.warp.vox = NaN;
                matlabbatch{5}.spm.spatial.preproc.warp.bb = [NaN NaN NaN
                    NaN NaN NaN];

                spm_jobman('run',matlabbatch)

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %Rename and save scan-specific segmentation results 
                destination_anat_c1 = replace(destination_anat,"sub-","c1sub-");
                destination_anat_c2 = replace(destination_anat,"sub-","c2sub-");
                destination_anat_c3 = replace(destination_anat,"sub-","c3sub-");
                destination_anat_c4 = replace(destination_anat,"sub-","c4sub-");
                destination_anat_c5 = replace(destination_anat,"sub-","c5sub-");
                
                movefile(destination_anat_c1,replace(destination_anat_c1,"Day1",rest_str));
                movefile(destination_anat_c2,replace(destination_anat_c2,"Day1",rest_str));
                movefile(destination_anat_c3,replace(destination_anat_c3,"Day1",rest_str));
                movefile(destination_anat_c4,replace(destination_anat_c4,"Day1",rest_str));
                movefile(destination_anat_c5,replace(destination_anat_c5,"Day1",rest_str));
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                copyfile(destination_func,[realignment_directory_out '/sub-S' num2str(sub) '_ses-Day' num2str(day) '_task-' rest_str '_acq-mb3me4_echo-' num2str(echo) '_bold.nii']);

                V_e1=spm_vol([realignment_directory_out '/sub-S' num2str(sub) '_ses-Day' num2str(day) '_task-' rest_str '_acq-mb3me4_echo-' num2str(echo) '_bold_realigned.nii']);

            elseif echo > 1

                copyfile(destination_func,[realignment_directory_out '/sub-S' num2str(sub) '_ses-Day' num2str(day) '_task-' rest_str '_acq-mb3me4_echo-' num2str(echo) '_bold_realigned.nii']);

                V=spm_vol(destination_func);
                Y=spm_read_vols(V);

                for i = 1:numel(V)
                    V(i).mat=V_e1(i).mat;
                    V(i).fname = [realignment_directory_out '/sub-S' num2str(sub) '_ses-Day' num2str(day) '_task-' rest_str '_acq-mb3me4_echo-' num2str(echo) '_bold_realigned.nii'];
                end

                for t = 1:numel(V)
                    spm_write_vol(V(t),Y(:,:,:,t));
                end
            end

        end
        rmdir(realignment_directory_temp, 's');
    end
end

