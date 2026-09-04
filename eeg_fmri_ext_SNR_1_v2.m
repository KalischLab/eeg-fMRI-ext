%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% eeg_fmri_ext_SNR_1_v2.m
%
% written by B.Meyer, Neuroimaging Center Mainz, 2026
% corrected version of eeg_fmri_ext_SNR_1.m
%
% Motion-corrects each echo of a multi-echo BOLD run and produces a GM mask
% in that run's functional space, as input to eeg_fmri_ext_SNR_2_v2.sh.
%
% All four echoes come from the same excitation, so the head is in the same
% position for all of them. Motion is therefore estimated ONCE, on echo 1
% (shortest TE -> least dropout, highest SNR for registration), and the
% resulting per-volume transforms are applied to all four echoes, which are
% then each resliced. This keeps the four echoes strictly comparable.
%
% Changes vs. v1:
%   1. realign ESTIMATE (not estwrite) on echo 1, so the per-volume motion
%      transforms live in the split files' headers. v1 used which = [0 1]
%      ("Mean Image Only"), which writes no resliced images at all, while
%      matlabbatch{3} then tried to concatenate "Realigned Images".
%   2. The per-volume affines are copied from the ESTIMATE stage. v1 copied
%      them from the resliced+concatenated echo-1 file, where every volume
%      shares one identical affine, so no motion information was transferred.
%   3. Echoes 2-4 are actually resliced. v1 only rewrote their headers, so
%      their voxel data stayed motion-contaminated while echo 1 was corrected.
%   4. All four echoes now get identical treatment (same interpolation, same
%      FOV masking, same output datatype), so echo effects are not confounded
%      with processing differences. Output is float32, not int16, to avoid
%      quantising the interpolated data.
%   5. Non-steady-state volumes can be discarded before realignment (n_dummy).
%   6. No state leaks between iterations. In v1, V_e1 was only assigned when
%      echo 1 succeeded and was never cleared, so a missing echo-1 or T1 file
%      caused the PREVIOUS run's/subject's headers to be applied silently.
%   7. The T1 is segmented once per subject and the tissue maps are carried
%      into each run's space via coreg "other", instead of re-segmenting the
%      same T1 for every run (6x less work, identical result).
%   8. TR is read from the BIDS sidecar and written into the 4D header, so
%      script 2 can detrend. v1 set cat.RT = NaN.
%   9. Realignment parameters are kept, for motion-vs-tSNR QC.
%
% Requires SPM12 (spm_file_merge with an RT argument, i.e. r6906 or later)
% and MATLAB R2016b or later (local functions in a script file).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all

spm('defaults','FMRI');
spm_jobman('initcfg');

%% ----------------------------- configuration -----------------------------

bids_dir = '/mnt/nicshare/EEG_fmri_ext/data/bids/Nifti/';
dir_temp = '/home/mey1b/EEG_fmri_ext/data/temp/';          % scratch, wiped per run
dir_seg  = '/home/mey1b/EEG_fmri_ext/data/seg_cache/';     % one segmentation per subject
dir_out  = '/home/mey1b/EEG_fmri_ext/data/EEG_fMRI_ext_realigned_for_tSNR/';
tpm_file = '/usr/local/MATLAB/Toolboxes/SPM12/spm/tpm/TPM.nii';

subjects  = 10:34;
restscans = 1:6;                 % 1, 6: without EEG;  2, 3, 4, 5: with EEG
rest_day  = [1 2 2 2 2 3];       % rest scan -> acquisition day
n_echoes  = 4;

% ---------------------------------------------------------------------------
% SET THIS. Number of non-steady-state ("dummy") volumes to discard before
% realignment. Volumes acquired before T1 equilibrium are much brighter than
% the rest of the series and inflate the temporal SD, which depresses tSNR.
% If your reconstruction already discards them, leave this at 0. Check the
% first few volumes of one run before trusting any absolute tSNR value.
n_dummy = 0;
% ---------------------------------------------------------------------------

overwrite = false;               % false: skip runs whose outputs already exist

%% ------------------------------- setup -----------------------------------

for d = {dir_out, dir_seg}
    if ~exist(d{1},'dir'), mkdir(d{1}); end
end

log_file = fullfile(dir_out, sprintf('realign_log_%s.txt', datestr(now,'yyyymmdd_HHMMSS')));
fid_log  = fopen(log_file,'w');
if fid_log < 0
    warning('cannot write log to %s; logging to console only', log_file);
    fid_log = 1;
end
if fid_log == 1
    logmsg = @(varargin) fprintf(1, varargin{:});
else
    logmsg = @(varargin) cellfun(@(f) fprintf(f, varargin{:}), {1, fid_log});
end

logmsg('n_dummy = %d\n\n', n_dummy);

%% ------------------------------ main loop --------------------------------

for sub = subjects

    sub_str   = sprintf('sub-S%d', sub);
    path_anat = fullfile(bids_dir, sub_str, 'ses-Day1', 'anat');
    anat_gz   = fullfile(path_anat, sprintf('%s_ses-Day1_T1w.nii.gz', sub_str));

    if exist(anat_gz,'file') ~= 2
        logmsg('SKIP %s: no T1 (%s)\n', sub_str, anat_gz);
        continue
    end

    % --- segment the T1 once per subject, in its native space -------------
    seg_dir  = fullfile(dir_seg, sub_str);
    anat_seg = fullfile(seg_dir, sprintf('%s_ses-Day1_T1w.nii', sub_str));
    c_seg    = arrayfun(@(c) fullfile(seg_dir, sprintf('c%d%s_ses-Day1_T1w.nii', c, sub_str)), ...
                        1:5, 'UniformOutput', false);

    if ~all(cellfun(@(f) exist(f,'file')==2, c_seg))
        if ~exist(seg_dir,'dir'), mkdir(seg_dir); end
        copyfile(anat_gz, seg_dir);
        gunzip(fullfile(seg_dir, sprintf('%s_ses-Day1_T1w.nii.gz', sub_str)));
        delete(fullfile(seg_dir, sprintf('%s_ses-Day1_T1w.nii.gz', sub_str)));

        logmsg('%s: segmenting T1\n', sub_str);
        try
            spm_jobman('run', segment_batch(anat_seg, tpm_file));
        catch ME
            logmsg('FAIL %s: segmentation error: %s\n', sub_str, ME.message);
            continue
        end
    end

    for rest = restscans

        day      = rest_day(rest);
        rest_str = sprintf('rest%d', rest);
        tag      = sprintf('%s_ses-Day%d_task-%s', sub_str, day, rest_str);

        path_func = fullfile(bids_dir, sub_str, sprintf('ses-Day%d', day), 'func');

        func_gz  = cell(n_echoes,1);
        func_nii = cell(n_echoes,1);
        f_out    = cell(n_echoes,1);
        for e = 1:n_echoes
            func_gz{e}  = fullfile(path_func, sprintf('%s_acq-mb3me4_echo-%d_bold.nii.gz', tag, e));
            func_nii{e} = sprintf('%s_acq-mb3me4_echo-%d_bold.nii', tag, e);
            f_out{e}    = fullfile(dir_out, sprintf('%s_acq-mb3me4_echo-%d_bold_realigned.nii', tag, e));
        end

        % One motion estimate is shared by all echoes, so a run is only
        % usable if every echo is present. Partial runs are skipped whole.
        if ~all(cellfun(@(f) exist(f,'file')==2, func_gz))
            logmsg('SKIP %s: not all %d echoes present\n', tag, n_echoes);
            continue
        end

        % Resume check must include the GM mask, not just the BOLD outputs:
        % skipping a run whose mask is missing would leave script 2 unable to
        % process it.
        f_gm = fullfile(dir_out, sprintf('c1%s_T1w.nii', tag));
        if ~overwrite && all(cellfun(@(f) exist(f,'file')==2, f_out)) ...
                && exist(f_gm,'file') == 2
            logmsg('SKIP %s: outputs already exist\n', tag);
            continue
        end

        logmsg('---- %s ----\n', tag);

        try
            %% ---- stage in scratch and split 4D -> 3D --------------------
            if exist(dir_temp,'dir'), rmdir(dir_temp,'s'); end
            mkdir(dir_temp);

            split_files = cell(n_echoes,1);
            for e = 1:n_echoes
                copyfile(func_gz{e}, dir_temp);
                gz_local = fullfile(dir_temp, sprintf('%s_acq-mb3me4_echo-%d_bold.nii.gz', tag, e));
                gunzip(gz_local);
                delete(gz_local);

                % Realignment needs one affine per volume, so split to 3D.
                Vo = spm_file_split(fullfile(dir_temp, func_nii{e}), dir_temp);
                f  = {Vo.fname}';
                clear Vo                                   % release the mapping

                if n_dummy > 0
                    if numel(f) <= n_dummy
                        error('run has only %d volumes, n_dummy = %d', numel(f), n_dummy);
                    end
                    for k = 1:n_dummy
                        delete(f{k});
                    end
                    f = f(n_dummy+1:end);
                end

                split_files{e} = f;
                delete(fullfile(dir_temp, func_nii{e}));   % the 4D copy
            end

            n_vol = cellfun(@numel, split_files);
            if ~all(n_vol == n_vol(1))
                error('echoes have different volume counts: %s', mat2str(n_vol));
            end
            n_vol = n_vol(1);
            logmsg('  %d volumes per echo\n', n_vol);

            %% ---- estimate motion ONCE, on echo 1 -----------------------
            % ESTIMATE only: this writes the per-volume rigid-body transforms
            % into each split file's header without touching the voxel data.
            mb = [];
            mb{1}.spm.spatial.realign.estimate.data             = split_files(1);
            mb{1}.spm.spatial.realign.estimate.eoptions.quality = 0.9;
            mb{1}.spm.spatial.realign.estimate.eoptions.sep     = 4;
            mb{1}.spm.spatial.realign.estimate.eoptions.fwhm    = 5;
            mb{1}.spm.spatial.realign.estimate.eoptions.rtm     = 1;   % two-pass, register to mean
            mb{1}.spm.spatial.realign.estimate.eoptions.interp  = 2;
            mb{1}.spm.spatial.realign.estimate.eoptions.wrap    = [0 0 0];
            mb{1}.spm.spatial.realign.estimate.eoptions.weight  = '';
            spm_jobman('run', mb);

            % Keep the motion parameters: needed to check whether tSNR
            % differences between days are really motion differences.
            rp = dir(fullfile(dir_temp,'rp_*.txt'));
            if isempty(rp)
                logmsg('  WARNING: no rp_*.txt produced\n');
            else
                copyfile(fullfile(dir_temp, rp(1).name), ...
                         fullfile(dir_out, sprintf('rp_%s.txt', tag)));
            end

            %% ---- apply that estimate to every echo, then reslice -------
            TR = read_tr(path_func, tag);
            if TR == 0
                logmsg('  WARNING: no RepetitionTime in sidecar; script 2 will not detrend\n');
            end

            for e = 1:n_echoes

                if e > 1
                    % Copy the motion-corrected orientation volume by volume.
                    % This is the step v1 got wrong: it took the affine from
                    % the resliced echo-1 output, where all volumes share one
                    % affine and the motion is therefore no longer encoded.
                    for t = 1:n_vol
                        spm_get_space(split_files{e}{t}, spm_get_space(split_files{1}{t}));
                    end
                end

                % Now actually resample the voxel data. v1 never did this for
                % echoes 2-4, so their motion was never removed.
                mb = [];
                mb{1}.spm.spatial.realign.write.data             = split_files{e};
                mb{1}.spm.spatial.realign.write.roptions.which   = [2 1];  % all images + mean
                mb{1}.spm.spatial.realign.write.roptions.interp  = 4;
                mb{1}.spm.spatial.realign.write.roptions.wrap    = [0 0 0];
                mb{1}.spm.spatial.realign.write.roptions.mask    = 1;
                mb{1}.spm.spatial.realign.write.roptions.prefix  = 'r';
                spm_jobman('run', mb);

                % Back to 4D. float32 keeps the interpolated values intact;
                % v1's int16 quantised echo 1 but not echoes 2-4.
                r_files = spm_file(split_files{e}, 'prefix', 'r');
                if ~all(cellfun(@(f) exist(f,'file')==2, r_files))
                    error('reslicing produced no r* files for echo %d', e);
                end
                spm_file_merge(r_files, f_out{e}, 16, TR);
                logmsg('  echo %d -> %s\n', e, f_out{e});
            end

            %% ---- GM mask in this run's functional space ----------------
            % The tissue maps are moved with the T1 by coreg "other", so the
            % single per-subject segmentation is reused for every run.
            mean_e1 = spm_select('FPList', dir_temp, '^mean.*echo-1_bold.*\.nii$');
            if isempty(mean_e1)
                error('no mean image found for echo 1');
            end
            mean_e1 = strtrim(mean_e1(1,:));

            anat_run = fullfile(dir_temp, sprintf('%s_ses-Day1_T1w.nii', sub_str));
            copyfile(anat_seg, anat_run);
            c_run = cell(5,1);
            for c = 1:5
                c_run{c} = fullfile(dir_temp, sprintf('c%d%s_ses-Day1_T1w.nii', c, sub_str));
                copyfile(c_seg{c}, c_run{c});
            end

            mb = [];
            mb{1}.spm.spatial.coreg.estimate.ref                = {[mean_e1 ',1']};
            mb{1}.spm.spatial.coreg.estimate.source             = {[anat_run ',1']};
            mb{1}.spm.spatial.coreg.estimate.other              = c_run;
            mb{1}.spm.spatial.coreg.estimate.eoptions.cost_fun  = 'nmi';
            mb{1}.spm.spatial.coreg.estimate.eoptions.sep       = [4 2];
            mb{1}.spm.spatial.coreg.estimate.eoptions.tol       = ...
                [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
            mb{1}.spm.spatial.coreg.estimate.eoptions.fwhm      = [7 7];
            spm_jobman('run', mb);

            % Run-specific names, built explicitly rather than by string
            % substitution. Script 2 reconstructs exactly these.
            for c = 1:5
                movefile(c_run{c}, fullfile(dir_out, sprintf('c%d%s_T1w.nii', c, tag)));
            end
            copyfile(mean_e1, fullfile(dir_out, sprintf('%s_meanEPI_echo-1.nii', tag)));

            rmdir(dir_temp,'s');

        catch ME
            logmsg('FAIL %s: %s\n', tag, ME.message);
            for k = 1:numel(ME.stack)
                logmsg('      at %s line %d\n', ME.stack(k).name, ME.stack(k).line);
            end
            if exist(dir_temp,'dir'), rmdir(dir_temp,'s'); end
            continue
        end
    end
end

logmsg('\ndone. log: %s\n', log_file);
fclose(fid_log);


%% ------------------------------ subfunctions -----------------------------

function TR = read_tr(path_func, tag)
% Read RepetitionTime from the BIDS JSON sidecar (echo 1).
% Returns 0 (the NIfTI convention for "not set") if unavailable: writing NaN
% into pixdim4 would make the TR unparseable for script 2.
TR = 0;
f = fullfile(path_func, sprintf('%s_acq-mb3me4_echo-1_bold.json', tag));
if exist(f,'file') == 2
    try
        j = spm_jsonread(f);
        if isfield(j,'RepetitionTime') && isscalar(j.RepetitionTime) ...
                && isfinite(j.RepetitionTime) && j.RepetitionTime > 0
            TR = double(j.RepetitionTime);
        end
    catch
        % leave 0
    end
end
end


function mb = segment_batch(anat_file, tpm_file)
% Unified segmentation, native-space tissue maps only. Same parameters as v1.
ngaus = [1 1 2 3 4 2];
mb = [];
mb{1}.spm.spatial.preproc.channel.vols     = {[anat_file ',1']};
mb{1}.spm.spatial.preproc.channel.biasreg  = 0.001;
mb{1}.spm.spatial.preproc.channel.biasfwhm = 60;
mb{1}.spm.spatial.preproc.channel.write    = [0 0];
for t = 1:6
    mb{1}.spm.spatial.preproc.tissue(t).tpm    = {sprintf('%s,%d', tpm_file, t)};
    mb{1}.spm.spatial.preproc.tissue(t).ngaus  = ngaus(t);
    mb{1}.spm.spatial.preproc.tissue(t).native = [double(t <= 5) 0];
    mb{1}.spm.spatial.preproc.tissue(t).warped = [0 0];
end
mb{1}.spm.spatial.preproc.warp.mrf     = 1;
mb{1}.spm.spatial.preproc.warp.cleanup = 1;
mb{1}.spm.spatial.preproc.warp.reg     = [0 0.001 0.5 0.05 0.2];
mb{1}.spm.spatial.preproc.warp.affreg  = 'mni';
mb{1}.spm.spatial.preproc.warp.fwhm    = 0;
mb{1}.spm.spatial.preproc.warp.samp    = 3;
mb{1}.spm.spatial.preproc.warp.write   = [0 0];
mb{1}.spm.spatial.preproc.warp.vox     = NaN;
mb{1}.spm.spatial.preproc.warp.bb      = [NaN NaN NaN; NaN NaN NaN];
end
