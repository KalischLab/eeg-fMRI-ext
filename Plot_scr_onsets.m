%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting SCR - Spontaneous Responses (Onset Counts)
% by M. Gudina, Neuroimaging Center, Mainz 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; close all;

% Load the data
data = readtable('Onset_Counts.csv');

% Extract the Onset values and grouping variable
onset_values = data.Onsets;
time_groups = categorical(data.Time);

% Plotting the data
figure('Name', 'Spontaneous Responses (Onset Counts)', 'Position', [100, 100, 400, 550]);
hold on; 

%%%%%%%%%%% PLOT BOX PLOTS %%%%%%%%%%
% Plot using the grouping variable
boxplot(onset_values, time_groups, 'Colors', 'k', 'Symbol', '', 'Widths', 0.4);

% Draw the individual dots for each resting state
groups = categories(time_groups);
for i = 1:length(groups)
    idx = (time_groups == groups{i});
    plot(i * ones(sum(idx), 1), onset_values(idx), 'ok', 'MarkerFaceColor', 'k', 'MarkerSize', 5, 'HandleVisibility', 'off');
end

% % Set the box to a different color (e.g., blue) and the median line to another (e.g., red)
% set(findobj(gca, 'Tag', 'Box'), 'Color', 'b');
 set(findobj(gca, 'Tag', 'Median'), 'Color', 'r');

% Remove whiskers and caps
set(findobj(gca, 'Tag', 'Whisker'), 'Visible', 'off'); 
set(findobj(gca, 'Tag', 'Upper Adjacent Value'), 'Visible', 'off');
set(findobj(gca, 'Tag', 'Lower Adjacent Value'), 'Visible', 'off');

%%%%%%%%%%%% FORMATTING Graph %%%%%%%%%%%%
set(gca, 'FontSize', 20);  
title('Spontaneous Skin Conductance Responses');
ylabel('nSCR');
xlabel('rsfMRI Runs');
ylim([-5, 50]);
xticklabels({'R1.1', 'R2.1', 'R2.2', 'R2.3', 'R2.4', 'R3.1'});

hold off;