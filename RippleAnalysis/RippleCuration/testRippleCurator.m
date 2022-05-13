% TestRippleCurator

%% Import libaries
% Buzcode GitHub Repo
addpath(genpath("C:\Users\Linus Meienberg\Documents\Buzcode"));
addpath(genpath('C:\Users\Linus Meienberg\Documents\MultiAreaBundleProject\RippleAnalysis'))

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


%% Load Ripple Events from manual annotation
load("C:\Users\Linus Meienberg\Documents\MultiAreaBundleProject\RippleAnalysis\RippleDetection\rTBY33S7_DS0505_manual.mat")
ripplesManual.clusterID = ripple_classes;
% 0 : SPWR
% 1 : Ripple
% 2 : false Positive
ripplesManual.timestamps = ripple_timestamps;
ripplesManual.n = size(ripplesManual.timestamps,1);
%% Load Feature Table from manual annotation
load("D:\INI\SemesterArbeitBaran\RippleAnalysis\RippleDetection\220504_DatasetConstruction\rTBY33S7_DS0511_manual_features.mat")
featureTableManual = featureTable;
ripplesManual.centers = table2array(featureTableManual(:,'center Timepoint'));

%% Convert Ripple ID to categorical array
ripple_categories = categorical(ripplesManual.clusterID,[0 1 2],{'SPWR','R','fp'});


%% Spin up RippleCurator
c = RippleCurator(2000,lfp.lowpass,ripplesManual.timestamps,ripplesManual.centers,ripple_categories,featureTableManual);