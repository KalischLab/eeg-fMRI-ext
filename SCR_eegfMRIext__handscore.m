%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EEG_fMRI_ext - SCR handscore (A. Gerlicher, 2018)
% export the final data to CSV file
% Updated by M. Gudina, Neuroimaging Center, Mainz 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
addpath('/mnt/nicshare/EEG_fmri_ext/code/milkiyas/SCR_analysis/handscore gui')
readdir = '/mnt/nicshare/EEG_fmri_ext/data/SC/SC_data';
logdir  = '/mnt/nicshare/EEG_fmri_ext/data/SC/SCR_log';
%plotdir = 'C:\Users\elena\Desktop\NIC_2025\Rep_AA\S3_data\log\';

files = dir(fullfile(readdir, 'sub_20*_exp_day*.mat')); 
for i = 1:length(files)
    fname = files(i).name;
    parts = split(erase(fname, ".mat"), ["sub_", "_exp_day"]);
    subject = str2double(parts{2}); day = str2double(parts{3});
    matfile = fullfile(readdir, fname);
%%%%%%%%%%%%%%%%%%%%%%%%%% Check if already analyzed %%%%%%%%%%%%%%%%%%%%%
    expected_output = strcat('SCR_Day', num2str(day), '_S_', num2str(subject), '_filtered.mat');
    if exist(expected_output, 'file')
        fprintf('Skipping Subject %d Day %d - Already done!\n', subject, day);
        continue; 
    end

%matfile     = ([readdir filesep 'sub_' num2str(subject) '_exp_day' num2str(day) '.mat' ]);
data        = importdata(matfile);
scorewindow = 7000;
resp_crit   = 0.02;

% fin         = fopen([readdir filesep 'S' num2str(subject) '_RepDopa_Day' num2str(day) '.txt']);
fin         = fopen([logdir filesep num2str(subject) '_eeg_fmri_ext_Day' num2str(day) '.txt']);
if fin == -1
    fprintf('Notice: No log file for Subject %d Day %d. Will score without sorting.\n', subject, day);
    has_log = false;
else
    has_log = true;
log_data    = textscan(fin, '%d %d %d %f %f','Headerlines', 3);
trialtype   = cell2mat(log_data(1,2));

if day == 1
    trialtype   = trialtype(1:20);
elseif day == 2
    trialtype   = trialtype(1:30);
elseif day == 3
    trialtype   = trialtype(1:20);
end;
end;
    
SCR_channel      = 1;
stimulus_channel = 6; 

%%%%%%%%%%%%%%%%%%%%%% filter data to reduce MRI artifact: %%%%%%%%%%%%%%%%%%%%%%

scr     = data.data(:,SCR_channel);
cutoff  = 1;
sampler = 1000;
[b,a]   = butter(2,cutoff/(sampler/2));
scr2    = filter(b,a,scr);
scr     = scr2(250:length(scr2)); % scr2(250:length(scr2))why this?


%%%%%%%%%%%%%%%%%%%%%% prepare data and call Victor's handscore_gui function %%%%%%%%%%%%

pin2        = data.data(:,stimulus_channel);
codeindices = find(pin2);
pin2idx     = codeindices(find([1 diff(codeindices)'~=1]));   

if pin2idx(1)== 1
   pin2idx  = pin2idx(2:length(pin2idx));
end;

% pin2idx = pin2idx(1:20); % remove 1 trial for subject 1
scrs = zeros(1,length(pin2idx));
ridx = randperm(length(pin2idx));

i = 1;
go = 1;

    while go
        x = pin2idx(ridx(i));
        [marker, scresponse] = hand_score_gui(scr(x:x+scorewindow)', 900, 4000, 7000);
        if marker == 1
            scrs(1,ridx(i)) = scresponse;
        elseif marker == 0
            scrs(1,ridx(i)) = 0;
        elseif marker == -1
            i = max(0,i-2);
        end
        i=i+1;
        if i > length(pin2idx)
            go = 0;
        end
    end

       
for criterioni = 1:length(scrs)
    if  scrs(criterioni) < resp_crit
        scrs(criterioni) = 0;
    end;  
end;

%%%%%%%%%%%%%%%%%%%%%%% relate SCRs and trialtypes here: %%%%%%%%%%%%%%%%%%%%%%
CS1 =[];
CS2 =[];

if has_log
 for orderi =1:length(scrs)
     if trialtype(orderi) == 111 ||trialtype(orderi) == 11  ||trialtype(orderi) == 1
        CS1 = [CS1 scrs(orderi)];
     elseif trialtype(orderi) == 222 || trialtype(orderi) == 22 || trialtype(orderi) == 2
        CS2 = [CS2 scrs(orderi) ]; 
     end;
end; 

%%%%%%%%%%%%%%%%%%%%%%% Please note, how many valid trials there are! %%%%%%%%%%%%%%%%%%%%%%

CS1_valid   = nnz(CS1);
CS2_valid   = nnz(CS2);
Total_valid = CS1_valid + CS2_valid
else
    CS1 = NaN; % If there is no log_file;
    CS2 = NaN;
    Total_valid = nnz(scrs);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot individual SCRs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure(1);
% 
% if day == 1 
%    x    = [1:length(CS1)];
%    plot(CS1,'r-d', 'LineWidth', 1.5)
%    hold on
%    plot(CS2, 'b-s', 'LineWidth',1.5)
%    %axis([0 10.5 0 3.5])% 3.5 as standard. if necessary, change axis here!!!day 1!! axis(minx maxx miny maxy)
%    set(gca, 'XTick', [1:10])
% elseif day == 2
%    x    = [1:length(CS1)];
%    plot(CS1,'r-d', 'LineWidth', 1.5)
%    hold on
%    plot(CS2, 'b-s', 'LineWidth',1.5)
%    set(gca,'XTick',[1:15])
%    %axis([0 15.5 0 3.5]) % change axis here for day 2!! axis(minx maxx miny maxy)
% elseif day == 3
%    Renewal   = [1:10];
%    Recovery  = [11:20];  
%    plot(Renewal, CS1(1:10),'r-d', 'LineWidth', 1.5)
%    hold on
%    plot(Renewal, CS2(1:10),'b-s', 'LineWidth', 1.5)
% %    plot(Recovery, CS1(11:20),'r-d', 'LineWidth', 1.5)
% %    plot(Recovery, CS2(11:20),'b-s', 'LineWidth', 1.5)
%    %axis([0 20.5 0 3.5]) %axis(minx maxx miny maxy)
%    set(gca, 'XTick', [1:16])
% elseif day == 4
%    Renewal   = [1:4];
%    Recovery  = [5:8];  
%    plot(Renewal, CS1(1:4),'r-d', 'LineWidth', 1.5)
%    hold on
%    plot(Renewal, CS2(1:4),'b-s', 'LineWidth', 1.5)
%    plot(Recovery, CS1(5:8),'r-d', 'LineWidth', 1.5)
%    plot(Recovery, CS2(5:8),'b-s', 'LineWidth', 1.5)
%    %axis([0 8.5 0 3.5]) %axis(minx maxx miny maxy)
%    set(gca, 'XTick', [1:8])   
% end;
% 
% legend('CS+', 'CS-')
% legend boxoff
% box('off')
% xlabel('Trials')
% ylabel('Skin Conductance Response')
 t = strcat('SCR_Day', num2str(day), '_S_', num2str(subject),'_filtered'); 
% title(t);
% cd(plotdir)
% saveas(1,t, 'jpg')
% cd(matdir)
save(t, 'CS1', 'CS2', 'Total_valid','scrs')
%dlmwrite('SCR_EEG_fMRI_Ext.csv', [subject, day, mean(CS1), mean(CS2)], '-append');
csv_file = fullfile(pwd, 'SCR_EEG_fMRI_Ext.csv');
    trials = [10, 15, 10]; 
    
    % 1. Create or Load the Grid using LOW-LEVEL commands ONLY
    C = cell(2, 71); 
    if exist(csv_file, 'file')
        fclose('all'); 
        fid = fopen(csv_file, 'r');
        if fid == -1
            error('\n[!] Cannot read CSV. Please close the file or delete it and run the script again.\n');
        end
        
        row = 1;
        while ~feof(fid)
            line = fgetl(fid);
            if ischar(line)
                cols = strsplit(line, ';', 'CollapseDelimiters', false);
                for c = 1:min(length(cols), 71)
                    val = strtrim(cols{c});
                    numVal = str2double(val);
                    if ~isnan(numVal) && ~strcmpi(val, 'NaN')
                        C{row, c} = numVal;
                    elseif ~isempty(val)
                        C{row, c} = val;
                    end
                end
                row = row + 1;
            end
        end
        fclose(fid);
    else
        % Build fresh headers if file doesn't exist
        C(1, [1, 2, 22, 52]) = {'Subject', 'Conditioning', 'Extinction', 'Retrieval'};
        C{2,1} = 'Sub_ID'; col = 2;
        for d = 1:3, for cond = {'CS+','CS-'}, for t = 1:trials(d), C{2,col} = sprintf('%s_T%d', cond{1}, t); col = col+1; end; end; end
    end
    
    % Find Row and Inject Data safely
    r = [];
    for r_idx = 3:size(C,1)
        if isnumeric(C{r_idx,1}) && C{r_idx,1} == subject
            r = r_idx; break;
        end
    end
    if isempty(r), r = size(C,1) + 1; C{r,1} = subject; end
    
    c_start = 2 + sum(trials(1:day-1)) * 2; t_day = trials(day);
    pad = @(x) [num2cell(x(1:min(end, t_day))), num2cell(NaN(1, t_day - min(length(x), t_day)))];
    
    C(r, c_start : c_start + 2*t_day - 1) = [pad(CS1), pad(CS2)];
    
    fid = fopen(csv_file, 'w');
    if fid == -1
        error('\n[!] Cannot write to CSV. Please make sure the file is not open in Excel!\n');
    end
    
    for row_idx = 1:size(C, 1)
        for col_idx = 1:size(C, 2)
            val = C{row_idx, col_idx};
            if ischar(val) || isstring(val)
                fprintf(fid, '%s', val);
            elseif isnumeric(val) && isscalar(val) && ~isnan(val)
                fprintf(fid, '%g', val);
            end
            if col_idx < size(C, 2)
                fprintf(fid, ';'); 
            else
                fprintf(fid, '\n'); 
            end
        end
    end
    fclose(fid);
 end
% 
% close all
