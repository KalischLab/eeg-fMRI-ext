% =========================================================================
% eegfMRIext_sc_rest_mg.m
% Author: M. Gudina, Neuroimaging Center Mainz 17.06.2026
% Preprocess resting-state SC data
% =========================================================================

close all; clear variables; clc;

%%%%%%%%%%%%%% DIRECTORIES %%%%%%%%%%%%%%
sublist   = 10:34; 
indir     = '/mnt/nicshare/EEG_fmri_ext/data/SC/SC_data';
outdir    = '/mnt/nicshare/EEG_fmri_ext/code/milkiyas/SCR_analysis'; 

if ~exist(outdir,'dir')
    mkdir(outdir);
end
%%%%%%%%%%%%%%% EXTRACTION SETTINGS %%%%%%%%%%%%%%
minutes_to_extract = 8; % Change this to 10 if you want the full 10 minutes!
sampling_rate = 1000;   % Raw sampling rate
samples_to_extract = minutes_to_extract * 60 * sampling_rate;

log_IDs      = cell(0,1);
log_Comments = cell(0,1);
%%%%%%%%%%%%%%% MAIN LOOP %%%%%%%%%%%%%%%%
for sub = sublist
    searchPattern = fullfile(indir, sprintf('sub_%d_*rest*.mat', sub));
    subjectFiles  = dir(searchPattern);
    
    for f = 1:length(subjectFiles)
        original_filename = subjectFiles(f).name;
        infile   = fullfile(indir, original_filename);
        outfile  = fullfile(outdir, original_filename);

        load(infile, 'data');
        sc_raw = data(:, 1);

        %find in data(:,6) onset of rsfMRI session
        %extract 6-10 min of data
        %do all subsequent processing steps on the extracted time window
        
        %find in data(:,6) onset of rsfMRI session
        onset = find(data(:, 6) > 0, 1, 'first');
        if isempty(onset)
            fprintf('WARNING: No trigger for %s. Skipping...\n', original_filename);
            log_IDs{end+1} = original_filename; log_Comments{end+1} = 'No trigger'; continue;
        end
        
        % Extract: 5000 samples before onset, 495000 samples after (8m + 15s)
        sc_raw = data(max(1, onset - 5000) : min(onset + 494999, size(data, 1)), 1);
        
%         % Filter & Downsample
%         sc = ft_preproc_lowpassfilter(sc_raw', 1000, 5)';
%         sc = sc(1:100:end);
%         t  = (0:length(sc)-1)' / 10;
        
        %lowpass filter (5Hz) using fieldtrip function
        sc_filtered = ft_preproc_lowpassfilter(sc_raw', 1000, 5)';
        sc = sc_filtered(1:100:end);
        t  = linspace(0, length(sc)/10, length(sc))';
        
        % Visual Inspection
        figure('Name', 'Rest Data Inspection', 'Position', [200, 200, 800, 400]);
        plot(t, sc, 'k', 'LineWidth', 1.5);
        title(['Filtered Rest Data: ', strrep(original_filename, '_', '\_')]);
        xlabel('Time (Seconds)'); grid on;
        
        r = input('Data OK? (1: Yes/Save, 0: Invalid/Discard): ');
        close;
        if r == 1
            save(outfile, 'sc', 't');
            saved_files{end+1} = original_filename;
            
      else if r == 0
            reason = input('  -> enter REASON for exclusion: ', 's');
            log_IDs{end+1} = original_filename;
            log_Comments{end+1} = reason;
        end
    end
end

%%%%%%%%%%%%%%% EXPORT EXCEL LOGS %%%%%%%%%%%%%%
if ~isempty(log_IDs)
    T_excl = table(log_IDs(:), log_Comments(:), 'VariableNames', {'Original_Filename', 'Reason'});
    writetable(T_excl, fullfile(outdir, 'Excluded_Rest_Log.xlsx'));
end
%%%%%%%%%%%%%% EXPORT SEPARATE EXCEL MATRIX (R1-R6)%%%%%%%%%%%%%%
rest_status = NaN(length(sublist), 6);

for i = 1:length(sublist)
    files = dir(fullfile(indir, sprintf('sub_%d_*rest*.mat', sublist(i))));
    for r = 1:min(length(files), 6)
        fname = files(r).name;
        if any(strcmp(saved_files, fname)), rest_status(i, r) = 1;
        elseif any(strcmp(log_IDs, fname)),  rest_status(i, r) = 0; end
    end
end

% % Build Table and instantly format Subject column
% T = array2table(rest_status, 'VariableNames', {'R1','R2','R3','R4','R5','R6'});
% T.Subject = compose('S%d', sublist(:)); % Instantly creates S10, S11, etc.
% T = T(:, [7, 1:6]); % Move the Subject column to the very front
% writetable(T, fullfile(outdir, 'RestingState_Summary.xlsx'));
end
