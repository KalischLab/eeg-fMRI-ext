%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EEG_fMRI_ext - Plotting SCR
% by M. Gudina, Neuroimaging Center, Mainz 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; close all;

addpath(genpath('/mnt/nicshare/EEG_fmri_ext/code/plotSpread/plotSpread')); 

% Load the data
data = readmatrix('SCR_EEG_fMRI_Ext.csv');
n = size(data, 1); 
mean_data = mean(data(:, 2:end), 1, 'omitnan');
std_data  = std(data(:, 2:end), 0, 1, 'omitnan');
sem_data  = std_data / sqrt(n);
% Calculate the mean and sem for Day 1 & 2
CS_plus_mean  = [mean_data(1:10), mean_data(21:35)];
CS_minus_mean = [mean_data(11:20), mean_data(36:50)];
CS_plus_sem  = [sem_data(1:10), sem_data(21:35)];
CS_minus_sem = [sem_data(11:20), sem_data(36:50)];
% Calculate the mean and sem for Day 3
subj_D3_plus  = mean(data(:, 52:61), 2, 'omitnan');
subj_D3_minus = mean(data(:, 62:71), 2, 'omitnan');
mean_D3_plus = mean(subj_D3_plus, 'omitnan');
sem_D3_plus  = std(subj_D3_plus, 'omitnan') / sqrt(n);
mean_D3_minus = mean(subj_D3_minus, 'omitnan');
sem_D3_minus  = std(subj_D3_minus, 'omitnan') / sqrt(n);
trials = 1:length(CS_plus_mean); % x-axis vector for Trial 1 to 25

% Plotting the data on a single graph
figure('Name', 'EegfMRI SCR : Conditioning; Extinction; Retrieval', 'Position', [100, 100, 950, 550]);
hold on; 

%%%%%%%%%% PLOT CONTINUOUS LINES (Days 1 & 2) %%%%%%%%%%%%%
% Plot CS- using errorbar (Mean +/- SEM)
errorbar(trials, CS_minus_mean, CS_minus_sem, ...
    '-bo', 'LineWidth', 1.5, 'MarkerFaceColor', 'b', 'DisplayName', 'CS- ');
% Plot CS+ using errorbar (Mean +/- SEM)
errorbar(trials, CS_plus_mean, CS_plus_sem, ...
    '-ro', 'LineWidth', 1.5, 'MarkerFaceColor', 'r', 'DisplayName', 'CS+ ');

%%%%%%%%%%% PLOT DAY 3 SUMMARY (Classic Box Plots) %%%%%%%%%%
% CS+ Boxplot (No dots)
boxplot(subj_D3_plus, 'Positions', 28, 'Colors', 'r', 'Symbol', '', 'Widths', 0.8);

% CS- Boxplot (No dots)
boxplot(subj_D3_minus, 'Positions', 30, 'Colors', 'b', 'Symbol', '', 'Widths', 0.8);

set(findobj(gca, 'Tag', 'Whisker'), 'Visible', 'off'); %remove the top whiskers
set(findobj(gca, 'Tag', 'Upper Adjacent Value'), 'Visible', 'off');
set(findobj(gca, 'Tag', 'Lower Adjacent Value'), 'Visible', 'off');

%%%%%%%%%%%% FORMATTING Graph %%%%%%%%%%%%
set(gca, 'FontSize', 20); % Increase axis and general font size
title('EegfMRI SCR : Conditioning; Extinction; Retrieval');
ylabel('SCR Amplitude (\muS)');
xlabel({' ', ' '});
xlim([1, 32]); 
ylim([-inf, 1.2]); 
xticks([1 5 10 15 20 25 28 30]);% Custom X-axis labels to clearly mark the trials vs the Day 
xticklabels({'1', '5', '10', '15', '20', '25', ' CS+', ' CS-'});
xline(10.5, '--k', 'HandleVisibility', 'off');%vertical lines to separate the experimental phases
xline(26.5, '-k', 'HandleVisibility', 'off');

text(0.145, -0.08, 'Day 1', 'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 14, 'Clipping', 'off');
text(0.565, -0.08, 'Day 2', 'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 14, 'Clipping', 'off');
text(0.903, -0.08, 'Day 3', 'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 14, 'Clipping', 'off');
legend('Location', 'best');
hold off;''