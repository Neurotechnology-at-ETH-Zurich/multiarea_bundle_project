% Demonstration workflow manual ripple annotation

%% add Buzcode Repo to search path
addpath(genpath('C:\Users\Linus\Documents\Buzcode'))
lfp_timestamps = (1:size(lfp_filtered,1))/2000;
%% import downsampled lfp traces of ripple active region
dHPC_ripple_channels = [120,85,84,86,83,87,82,88];
%lfp = ImporterDAT_multi('D:\INI\SemesterArbeitBaran\rTBY33\7_freely_behav_220315_145519\amplifier_ds.dat',256,dHPC_ripple_channels);
%lfp = ImporterDAT_multi('/media/baran/Linus/SemesterArbeitBaran_copy/rTBY33/7_freely_behav_220315_145519/amplifier_ds.dat',256,dHPC_ripple_channels);
lfp = ImporterDAT_multi('D:/SemesterArbeitBaran_copy/rTBY33/7_freely_behav_220315_145519/amplifier_ds.dat',256,dHPC_ripple_channels);
%% filter traces
% design a digital band-pass filter for the frequency range from 130Hz to 245Hz
% This isolates the lfp events in the 
d=designfilt('bandpassfir','FilterOrder',600,'StopbandFrequency1',125,'PassbandFrequency1',130,'PassbandFrequency2',200,'StopbandFrequency2',205,'SampleRate',2000);
% freqz(d,600) % visualize filter transfer function 
% perform zero phase digital filtering => no phase distortions in the
% signal, end effect equals to a filter with a transfer function equal to
% the squared magnitude of the original filter
lfp_detrended = detrend(double(lfp'));
lfp_filtered = filtfilt(d,lfp_detrended);

%% Generate automatic ripple annotations
[ripples_auto] = bz_FindRipples(lfp_detrended(:,2),lfp_timestamps, 'durations',[20 100],'minDuration',5,'frequency',2000,'passband',[180 220],'EMGThresh',0); %Buzsaki uses 20ms for min duration

%% Use Hilbert transform on bandpass filtered trace 86
lfp_filtered_hilbert = abs(hilbert(lfp_filtered));
%% Invoke LFPViewer on generated annotations
ripple_classes = zeros(size(ripples_auto.timestamps,1),1);
viewer = LFPViewer([lfp_detrended lfp_filtered lfp_filtered_hilbert*4],2000,ripples_auto.timestamps,ripple_classes,'rTBY33S7_ripples_auto');

%% Load manual ripple annotations
load("rTBY33S7_ripples_manual.mat")
%% Invoke LFPViewer
viewer = LFPViewer([lfp_detrended lfp_filtered lfp_filtered_hilbert (lfp_filtered_hilbert>33)*50],2000,ripples,'rTBY33S7_ripples_manual');
%% Load manual ripple annotations (updated)
load("rTBY33S7_ripples_manual.mat")
%% Save back to project folder using Buszaki compatible format (event struct)
n_ripples = size(ripples,1);
ripples_manual.timestamps = ripples;
ripples_manual.peaks = mean(ripples,2); % just take center point for now
ripples_manual.amplitude = ones(n_ripples,1); % just assign unit amplitude to each event
ripples_manual.amplitudeUnits = 'uV';
ripples_manual.eventID = zeros(n_ripples,1); % assign all ripples to class 0
ripples_manual.detectorinfo = 'Ripples annotated manualy to rTBY33 Session 7';
%save('D:\INI\SemesterArbeitBaran\rTBY33\7_freely_behav_220315_145519\amplifier.ripplesManual.events.mat','ripples_manual');