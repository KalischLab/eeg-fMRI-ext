%EegfMRIext_sc_rest_analysis

%Author: Benjamin Meyer, Neuroimaging Center Mainz, 23.06.2020
%Updated by: Milkiyas Gudina, 22.06.2026

close all; clear all; clc;
%addpath(genpath('/home/intern1/Documents/MATLAB/ledalab-349/'));

target_folder = '/mnt/nicshare/EEG_fmri_ext/code/milkiyas/SCR_analysis/Ledalabrest2/';
Ledalab(target_folder,'open','mat','smooth',{'mean',5},'analyze','CDA', 'optimize',2, 'export_era',[0 480 .02 1],'export_scrlist', [.02    1])
