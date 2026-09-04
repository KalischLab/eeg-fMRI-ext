% =========================================================================
% Extract NSCRs (Number of Onsets) from Ledalab scrlist Files
% Author M. Gudina, 12.08.2026, NiC, Mainz
% =========================================================================
close all; clear variables; clc;

%%%%%%%%%%%%%%%%%%%%% Directories %%%%%%%%%%%%%%%%%%%%%%%%
input_dir  = '/mnt/nicshare/EEG_fmri_ext/code/milkiyas/SCR_analysis/Ledalabrest2'; 
output_csv = fullfile(input_dir, 'Onset_Counts.csv');

files = dir(fullfile(input_dir, '*scrlist.mat'));
if isempty(files)
    error('[!] No scrlist.mat files found. Please check your input directory.');
end

% --- MAPPING LOGIC ---
% Format: [Rest_Number, Day_Number, Target_R_Column]
map_logic = [
    1, 1, 1;
    1, 2, 2;
    2, 2, 3;
    3, 2, 4;
    4, 2, 5;
    1, 3, 6
];

Sub_IDs     = [];
Rest        = []; 
Onset_Count = [];

%%%%%%%%%%%%%%%% Automatic loop for every file %%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:length(files)
    filename = files(i).name;
    filepath = fullfile(input_dir, filename);
    
    tokens = regexp(filename, 'sub_(\d+).*?rest_?(\d+).*?day_?(\d+)', 'tokens');
    
    if ~isempty(tokens)
        sub_num  = str2double(tokens{1}{1});
        rest_num = str2double(tokens{1}{2});
        day_num  = str2double(tokens{1}{3});
        
        match = find(map_logic(:,1) == rest_num & map_logic(:,2) == day_num);
        if isempty(match)
            continue; 
        end
        r_idx = map_logic(match, 3);
    else
        continue;
    end
    
    % Safely extract the onset array length
    try
        fileData = load(filepath); 
        
        if isfield(fileData, 'scrlist')
            target_var = fileData.scrlist;
        elseif isfield(fileData, 'scrList')
            target_var = fileData.scrList;
        else
            target_var = [];
        end
        
        if isempty(target_var)
            num_onsets = NaN;
        else
            num_onsets = 0; % Default to 0 true onsets
            
            % Check TTP method first, fallback to CDA method if needed
            if isfield(target_var, 'TTP') && isfield(target_var.TTP, 'onset')
                num_onsets = length(target_var.TTP.onset);
            elseif isfield(target_var, 'CDA') && isfield(target_var.CDA, 'onset')
                num_onsets = length(target_var.CDA.onset);
            end
        end
        
    catch
        num_onsets = NaN;
    end
    
    Sub_IDs(end+1, 1)     = sub_num;
    Rest(end+1, 1)        = r_idx;
    Onset_Count(end+1, 1) = num_onsets;
end

%%%%%%%%%%%%%%%% Format Data into 3 Columns %%%%%%%%%%%%%%%%%%%%%%%%
% Assign values to appropriately named column variables
ID = Sub_IDs;
Onsets = Onset_Count;

% Convert the numeric Rest array (1-6) into string labels
label_map = ["R1.1"; "R2.1"; "R2.2"; "R2.3"; "R2.4"; "R3.1"];
Time = label_map(Rest);

% Compile into the final 3-column table
ResultsTable = table(ID, Onsets, Time);

% Sort the table numerically by ID, then by Time to keep subjects grouped cleanly
ResultsTable = sortrows(ResultsTable, {'ID', 'Time'});

writetable(ResultsTable, output_csv);
fprintf('\nProcessing complete. Data saved to:\n%s\n', output_csv);