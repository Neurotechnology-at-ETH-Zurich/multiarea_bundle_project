%% DatasetConstruction
%
% This script demonstrates how to detect ripple events in the local field
% potential data and how the detected events can be inspected and modified
% using the LFPViewer
%
% May 2022, Linus Meienberg

%% Import Dependencies
addpath(genpath('C:\Users\Linus Meienberg\Documents\MultiAreaBundleProject\RippleAnalysis'))
%% Import downsampled electrophysiological recording
% for analysis of local field potential events, the high frequency content
% of the electrophysiological recordings can be ignored. We work with a
% downsampled recording file which has a samplerate of 2kHz
dHPC_ripple_channels = [120,85,84,86,83,87,82,88];
lfp = ImporterDAT_multi('D:\INI\SemesterArbeitBaran\rTBY33\7_freely_behav_220315_145519\amplifier_ds.dat',256,dHPC_ripple_channels);
%lfp = ImporterDAT_multi('/media/baran/Linus/SemesterArbeitBaran_copy/rTBY33/7_freely_behav_220315_145519/amplifier_ds.dat',256,dHPC_ripple_channels);
%lfp = ImporterDAT_multi('F:/SemesterArbeitBaran_copy/rTBY33/7_freely_behav_220315_145519/amplifier_ds.dat',256,dHPC_ripple_channels);
%% filter traces
% design a digital band-pass filter for the frequency range from 130Hz to 245Hz
% This isolates the lfp events in the 
d=designfilt('bandpassfir','FilterOrder',600,'StopbandFrequency1',125,'PassbandFrequency1',130,'PassbandFrequency2',245,'StopbandFrequency2',250,'SampleRate',2000);
% freqz(d,600) % visualize filter transfer function 
% perform zero phase digital filtering => no phase distortions in the
% signal, end effect equals to a filter with a transfer function equal to
% the squared magnitude of the original filter
lfp_detrended = detrend(double(lfp'));
lfp_filtered = filtfilt(d,lfp_detrended);

%% Generate automatic ripple annotations
% add Buzcode Repo to search path
addpath(genpath('C:\\Users\\Linus Meienberg\\Documents\\Buzcode'));
lfp_timestamps = (1:size(lfp_filtered,1))/2000;
%% Detect Ripples using different detector settings
[ripples_auto_hi] = bz_FindRipples(lfp_detrended(:,1),lfp_timestamps, 'durations',[20 100],'minDuration',5,'frequency',2000,'passband',[180 220],'EMGThresh',0); %Buzsaki uses 20ms for min duration
[ripples_auto_low] = bz_FindRipples(lfp_detrended(:,1),lfp_timestamps, 'durations',[20 100],'minDuration',5,'frequency',2000,'passband',[120 220],'EMGThresh',0);

%% Detect Ripples from bandpass power signal
lfp_filtered_hilbert = abs(hilbert(lfp_filtered));
lfp_filtered_hilbert_mean = mean(lfp_filtered_hilbert,2);
% histogram(lfp_filtered_hilbert_mean,100) % inspect bandpass power levels
detection_meanPowerPeakThreshold = 40;
detection_minPeakSeparationFrames = 0.01 * 2000;
detection_minPeakWidth = 0.010 * 2000; % minimum peak width / event duration
[meanPowerPeaksAmplitude, meanPowerPeaksLocation, meanPowerPeakWidth, ~] = findpeaks(lfp_filtered_hilbert_mean,'MinPeakDistance',...
    detection_minPeakSeparationFrames,'MinPeakHeight',detection_meanPowerPeakThreshold,...
    'MinPeakWidth',detection_minPeakWidth);
% calculate candidate intervals using the detected peak widths
n_candidates = numel(meanPowerPeaksAmplitude);
duration = min(meanPowerPeakWidth,150); % restrict half width to 300 frames (150ms@2kHz) max
candidate_intervals = [ meanPowerPeaksLocation - duration, meanPowerPeaksLocation + duration ] / 2000;
% merge overlapping and close events
detection_eventMergeTolerance = 0.01; % merge events that begin up to tolerance ms after prior event ended

merge_candidates = candidate_intervals(2:end,1) - candidate_intervals(1:end-1,2) < detection_eventMergeTolerance;

ripples_peakPower = [];
current_interval = candidate_intervals(1,:);
for i = 2:n_candidates
    if merge_candidates(i-1)
        current_interval = [current_interval(1) candidate_intervals(i,2)];
    else
        ripples_peakPower = [ripples_peakPower ; current_interval];
        current_interval = candidate_intervals(i,:);
    end
end
ripples_peakPower = [ripples_peakPower ; current_interval];

%% Delete intersecting events with priority ripples_auto_hi > ripples_auto_low > ripples_peakPower
% Delete intervals that intersect with ripples_auto_hi
for i = 1:size(ripples_auto_hi.timestamps,1)
    interval = ripples_auto_hi.timestamps(i,:);
    % look for overlap to ripples_auto_low
    % intervals without intersection end before the reference intervals
    % start or start after the reference interval ends
    to_delete = ~(ripples_auto_low.timestamps(:,2) < interval(1) | interval(2) < ripples_auto_low.timestamps(:,1));
    ripples_auto_low.peakNormedPower(to_delete) = [];
    ripples_auto_low.peaks(to_delete) = [];
    ripples_auto_low.timestamps(to_delete,:) = [];
    % look for overlap to ripples_peak_power
    to_delete = ~(ripples_peakPower(:,2) < interval(1) | interval(2) < ripples_peakPower(:,1));
    ripples_peakPower(to_delete,:) = [];
end

% Delete intervals that intersect with ripples_auto_low
for i = 1:size(ripples_auto_low.timestamps,1)
    interval = ripples_auto_low.timestamps(i,:);
    % look for overlap to ripples_peak_power
    to_delete = ~(ripples_peakPower(:,2) < interval(1) | interval(2) < ripples_peakPower(:,1));
    ripples_peakPower(to_delete,:) = [];
end

%% Assemble detected Ripples to dataset
n_ripples_tot = size(ripples_auto_hi.timestamps,1) + size(ripples_auto_low.timestamps,1) + size(ripples_peakPower,1);
ripple_timestamps = [ripples_auto_hi.timestamps ; ripples_auto_low.timestamps; ripples_peakPower];
ripple_classes = [zeros(size(ripples_auto_hi.timestamps,1),1) ; ones(size(ripples_auto_low.timestamps,1),1) ; 2*ones(size(ripples_peakPower,1),1) ];
ripple_classes = categorical(ripple_classes,[0 1 2],{'bzHigh','bzLow','pp'});

%% Invoke LFPViewer on generated annotations
centers = sum(ripple_timestamps,2)/2;
viewer = LFPViewer([lfp_detrended lfp_filtered lfp_filtered_hilbert*4],2000,ripple_timestamps,centers,ripple_classes,'rTBY33S7_DS0505_auto',false);
%% Invoke LFPViewer for manual annotation
viewer = LFPViewer([lfp_detrended lfp_filtered lfp_filtered_hilbert*4],2000,[],[],'rTBY33S7_DS0505_manual');
