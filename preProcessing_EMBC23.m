%-------------------------------------------------------------------------%
% Script: preProcessing_EMBC23.m
% Author: Remi Dagenais
% Date:   07/02/2023
% Descr:  Preprocess ABP and PPG signal from the calibration files.
%         Requires a peak finder detection(e.g. PhysioPeaksFinder - 
%         https://github.com/DagenaisR/PhysioPeaksFinder)
%         Run one section at a time and follow the instructions. 
% Ref:    Dagenais R., Mitsis G.D., Non-invasive estimation of arterial 
%         blood pressure fluctuations using a peripheral photoplethysmograph
%         inside the MRI scanner. EMBC23
%-------------------------------------------------------------------------%
clear variables; clc; clf;

% Identification
sub = 0000;
ses = 1;
sequence = 'rest';

%% Case by case processing (Load data, crop/sync data, resample etc.)
% Import PPG and ABP data (fill this part)
% sync data
% resample ABP and PPG at same Fs

%% Filtering PPG
% sampling rate
Fs = 250; %change this to actual sampling rate

% Filter PPG Signal into LP and HP [0.5-10 Hz]
[HPnum,HPdenum] = butter(4,0.5/Fs,'high');
[LPnum,LPdenum] = butter(4,10/Fs);
ppg_tmp = filtfilt(LPnum,LPdenum,ppg); ppg_filt = filtfilt(HPnum,HPdenum,ppg_tmp);

%% Detect PPG P1 (systolic peaks) and P3 (diastolic peaks)
PhysioPeaksFinder %Or any other peak detection algorithm

%% Run after exporting PhysioPeaksFinder PPG P1 & P3 idx
p3p1_idx = idx; clear idx;
p1idx    = p3p1_idx(1:2:end-1);
p3idx    = p3p1_idx(2:2:end);

%% Detect P1 (systolic peaks) from ABP signal
PhysioPeaksFinder

%% Run after exporting PhysioPeaksFinder ABP P1 idx
pp_idx = idx; clear idx;

%% Save pre-processed variables
filename = sprintf('SUB%i_SES0%i_%s',sub,ses,sequence);
save([filename '_processed.mat'],'sub','ses','sequence','time','abp','ppg_filt','pp_idx','p1_idx','p3_idx','Fs');

