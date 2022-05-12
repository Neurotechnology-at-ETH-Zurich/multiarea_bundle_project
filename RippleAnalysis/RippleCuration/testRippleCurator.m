% TestRippleCurator

%% Import libaries
% Buzcode GitHub Repo
addpath(genpath("C:\Users\Linus Meienberg\Documents\Buzcode"));
addpath(genpath('D:\INI\SemesterArbeitBaran\RippleAnalysis'))

%% Script Parameters
output_dir = "D:\INI\SemesterArbeitBaran\rTBY33\7_freely_behav_220315_145519";
lfp.filename = 'D:\INI\SemesterArbeitBaran\rTBY33\7_freely_behav_220315_145519\amplifier_ds.dat';
lfp.samplerate = 2000;
lfp.num_channels = 256;
lfp.ripple_channels = [121,120,85,84,86,83,87,82,88,119]; % channel number of contigous sites that show ripples
lfp.main_channel = 1; % index of the main channel for ripple detection in lfp.ripple_channels

%% Import LFP data
lfp.data = detrend(double(ImporterDAT_multi(lfp.filename,lfp.num_channels,lfp.ripple_channels))'); % detrend and subsequent op require (time,channels) format
lfp.timestamps = (1:size(lfp.data,1))/lfp.samplerate;

%% Filter Data
bandpass = designfilt('bandpassfir','FilterOrder',600,'StopbandFrequency1',125,'PassbandFrequency1',130,'PassbandFrequency2',200,'StopbandFrequency2',205,'SampleRate',2000);
lowpass = designfilt('lowpassfir','FilterOrder',600,'PassbandFrequency',245,'StopbandFrequency',250,'SampleRate',2000);
lfp.bandpass = filtfilt(bandpass,lfp.data);
lfp.lowpass = filtfilt(lowpass,lfp.data);


%% Load Ripple Events from automatic annotation
load("D:\INI\SemesterArbeitBaran\RippleAnalysis\RippleDetection\220504_DatasetConstruction\rTBY33S7_DS0505_auto_v2.mat")
ripples.clusterID = ripple_classes+1;
% 1 : bz_findRipples high
% 2 : bz_findRipples low
% 3 : peakPowerMethod
ripples.timestamps = ripple_timestamps;
ripples.n = size(ripples.timestamps,1);
%% Load Feature Table from automatic annotation
load("D:\INI\SemesterArbeitBaran\RippleAnalysis\RippleDetection\220504_DatasetConstruction\rTBY33S7_DS0511_auto_v2_features.mat")
featureTableAuto = featureTable;
ripples.centers = table2array(featureTable(:,'center Timepoint'));
%% Load Ripple Events from manual annotation
load("D:\INI\SemesterArbeitBaran\RippleAnalysis\RippleDetection\220504_DatasetConstruction\rTBY33S7_DS0505_manual.mat")
ripplesManual.clusterID = ripple_classes+1;
% 4 : SPWR
% 5 : Ripple
% 6 : false Positive
ripplesManual.timestamps = ripple_timestamps;
ripplesManual.n = size(ripplesManual.timestamps,1);
%% Load Feature Table from manual annotation
load("D:\INI\SemesterArbeitBaran\RippleAnalysis\RippleDetection\220504_DatasetConstruction\rTBY33S7_DS0511_manual_features.mat")
featureTableManual = featureTable;
ripplesManual.centers = table2array(featureTableManual(:,'center Timepoint'));

%%


%% Spin up RippleCurator
c = RippleCurator(2000,lfp.lowpass,ripplesManual.timestamps,ripplesManual.centers,ripplesManual.clusterID,featureTableManual);