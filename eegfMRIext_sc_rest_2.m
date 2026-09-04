% =========================================================================
% eegfMRIext_sc_rest_AutoRescue.m
% Reads your approved log, processes the raw data, and exports to Ledalab
% =========================================================================
close all; clear variables; clc;

%%%%%%%%%%%%% DIRECTORIES %%%%%%%%%%%%%%%%%%%
raw_dir   = '/mnt/nicshare/EEG_fmri_ext/data/SC/SC_data'; 
step1_dir = '/mnt/nicshare/EEG_fmri_ext/code/milkiyas/SCR_analysis';
outdir    = fullfile(step1_dir, 'Ledalabrest2');

if ~exist(outdir,'dir')
    mkdir(outdir);
end

files = dir(fullfile(raw_dir,'sub_*_rest*.mat'));
%%%%%%%%%%%%%% AUTO-PROCESS LOOP %%%%%%%%%%%%%%%%%%%%%
for i = 1:length(files)
    orig_name = files(i).name;
    infile = fullfile(raw_dir, orig_name);
    
    if exist(infile, 'file')
        % Load the RAW data matrix
        load(infile, 'data'); 
        
        % Find Onset & Extract (8m + buffers)
        onset = find(data(:, 6) > 0, 1, 'first');
        if isempty(onset)
            fprintf('Skipping %s (No trigger found)\n', orig_name); continue;
        end
        sc_raw = data(max(1, onset - 5000) : min(onset + 494999, size(data, 1)), 1);
        
        % Filter & Downsample to 10 Hz
        sc_filtered = ft_preproc_lowpassfilter(sc_raw', 1000, 5)';
        sc = sc_filtered(1:100:end);
        
        % Format for Ledalab
        perfect_time = (0:length(sc)-1) / 10;
          
        leda_data.conductance = sc(:)'; 
        leda_data.time        = perfect_time(:)';  
        leda_data.timeoff     = 0;
        leda_data.event = []; 
        leda_data.event(1).time = 0;
        leda_data.event(1).nid = 1;
        leda_data.event(1).name = 'rest_start';
        leda_data.event(1).userdata = [];
     
        % Overwrite the 'data' variable with Ledalab structure and save
        data = leda_data; 
        outfile = fullfile(outdir, ['ledalab_' orig_name]);
        save(outfile, 'data', '-v7');
        
        fprintf('Successfully processed and formatted: %s\n', ['ledalab_' orig_name]);
    else
        fprintf('WARNING: Could not find raw file %s\n', orig_name);
    end
end

fprintf('\n All files are ready for Ledalab.\n');