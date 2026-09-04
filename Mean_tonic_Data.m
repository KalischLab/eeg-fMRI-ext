%=========================================================================
% Extract Mean Tonic SCL and Format into R1-R6 Matrix
% Author M. Gudina, NIC, Mainz 12.08.2026, 
% =========================================================================
close all; clear variables; clc;
%%%%%%%%%%%%%%%%%%%%% Directories %%%%%%%%%%%%%%%%%%%%%%%%
input_dir = '/mnt/nicshare/EEG_fmri_ext/code/milkiyas/SCR_analysis/Ledalabrest2'; 
output_csv = fullfile(input_dir, 'Mean_Tonic_Data.csv');
files = dir(fullfile(input_dir, 'ledalab_sub_*_rest*.mat'));
% Mapping logic Format: [Rest_Num, Day_Num, Target_R_Column]
map_logic = [
    1, 1, 1;  % Rest 1, Day 1 -> R1
    1, 2, 2;  % Rest 1, Day 2 -> R2
    2, 2, 3;  % Rest 2, Day 2 -> R3
    3, 2, 4;  % Rest 3, Day 2 -> R4
    4, 2, 5;  % Rest 4, Day 2 -> R5
    1, 3, 6   % Rest 1, Day 3 -> R6
];
% Initialize empty arrays
Sub_IDs    = [];
Rest       = []; 
Mean_Tonic = [];
%%%%%%%%%%%%%%%% Automatic loop for every file %%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:length(files)
    filename = files(i).name;
    
    if contains(filename, '_era.mat') || contains(filename, '_scrlist.mat')
        continue;
    end
    
    filepath = fullfile(input_dir, filename);
    tokens = regexp(filename, 'sub_(\d+).*?rest_?(\d+).*?day_?(\d+)', 'tokens');
    
    if ~isempty(tokens)
        sub_num  = str2double(tokens{1}{1});
        rest_num = str2double(tokens{1}{2});
        day_num  = str2double(tokens{1}{3});
        
        % Check the lookup matrix for a match
        match = find(map_logic(:,1) == rest_num & map_logic(:,2) == day_num);
        if isempty(match)
            fprintf('[!] Unmapped combination (Rest %d, Day %d) in %s\n', rest_num, day_num, filename);
            continue; 
        end
        r_idx = map_logic(match, 3);
        
    else
        fprintf('[!] Filename format missing Day or Rest: %s\n', filename);
        continue;
    end
    
    % Load the 'analysis' struct from the .mat file
    try
        load(filepath, 'analysis'); 
        if exist('analysis', 'var') && isfield(analysis, 'tonicData')
            mean_val = mean(analysis.tonicData, 'omitnan');
        else
            fprintf('  -> Warning: No tonicData found in %s\n', filename);
            mean_val = NaN;
        end
    catch
        fprintf('  -> Error loading file: %s\n', filename);
        mean_val = NaN;
    end
    
    Sub_IDs(end+1, 1)    = sub_num;
    Rest(end+1, 1)       = r_idx;
    Mean_Tonic(end+1, 1) = mean_val;
    
    clear analysis; 
end

%%%%%%%%%%%%%%%% Format Data into 3 Columns %%%%%%%%%%%%%%%%%%%%%%%%
% Assign values to appropriately named column variables
ID = Sub_IDs;
SCL = Mean_Tonic;

% Convert the numeric Rest array (1-6) into string labels
label_map = ["R1.1"; "R2.1"; "R2.2"; "R2.3"; "R2.4"; "R3.1"];
Time = label_map(Rest);

% Compile into the final 3-column table
ResultsTable = table(ID, SCL, Time);

% Sort the table numerically by ID, then by Time to keep subjects grouped cleanly
ResultsTable = sortrows(ResultsTable, {'ID', 'Time'});

% Write to CSV
writetable(ResultsTable, output_csv);
fprintf('\n 3-Column data exported to:\n%s\n', output_csv);