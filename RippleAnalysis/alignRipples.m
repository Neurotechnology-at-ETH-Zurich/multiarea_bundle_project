function alignedCenterFrames  = alignRipples(signal,rippleCenterSeeds)
% alignRipples refine alignment of ripples based on cross correlation with
% template.
%
% signal            (:,1)       signal in which ripple events were
% detected. Assumed to be sampled at 2kHz.
% rippleCenterSeeds (ripples,1) list of estimated ripple center frames.
% Ripples are assumed to be 200Hz oscillations with a bell shaped
% amplitude.
%
% alignedCenterFrames (ripples,) the center frames of the aligned
% eventWindows

% start by using the provided seeds as central frames
rippleEpisodes.centralFrames = rippleCenterSeeds;

% extract centered windows around the seeds for each ripple
eventWindowFrames = 301;
eventWindowHalfFrames = floor(eventWindowFrames/2);
num_ripples = size(rippleCenterSeeds,1);
for i = 1:num_ripples
    startFrame = rippleEpisodes.centralFrames(i) - eventWindowHalfFrames;
    endFrame = rippleEpisodes.centralFrames(i) + eventWindowHalfFrames;
    rippleEpisodes.episodes(i,:) = signal(startFrame:endFrame);
end

% Refine the center frames by maximizing cross correlation to a
% signal template.
% Restrict the shift to +/- half a ripple cycle
% ripples ~= 200 Hz  = 10 frame cycle at 2000Hz sampling rate
maxlag = 10; % < half cycle of ripple signal in frames for stability
gamma = 3; % exponent controlls steepness of decay as shift -> maxlag

% Construct the alignment template using a bell shaped curve centered at
% the middle of the alignment window with a decay length of 30 frames (15ms
% @ 2kHz) multiplied with a sinusoidal oscillation at ~200Hz (~10
% frames/cycle @ 2kHz)
window = normalize( normpdf(-eventWindowHalfFrames:eventWindowHalfFrames,0,30),'range')';
template = sin(2*pi/12*(1:eventWindowFrames))' .* window;

% greedily align signals to template by shifting them up to maxlag
for i = 1:num_ripples
    % find signal offset by maximizing cross correlation
    [r,lags] = xcorr(template,rippleEpisodes.episodes(i,:),maxlag);
    % tamper cross correlation where lag -> +/- maxlag using a power
    % function
    corr_tamper = 1-(abs(lags')/maxlag).^gamma;
    [maxr,d] = max(r.*corr_tamper); % get highest cross correlation and its associated delay
    delay = d - maxlag - 1; % negative d => pad ripple A, positive d => pad rippleB
    offset = -delay;

    % shift the central frame of the aligned event
    rippleEpisodes.centralFrames(i) = rippleEpisodes.centralFrames(i) + offset;
    % extract the updated episode if the episode was shifted
    if offset ~= 0
        startFrame = rippleEpisodes.centralFrames(i) - eventWindowHalfFrames;
        endFrame = rippleEpisodes.centralFrames(i) + eventWindowHalfFrames;
        rippleEpisodes.episodes(i,:) = signal(startFrame:endFrame);
    end
end
    
alignedCenterFrames = rippleEpisodes.centralFrames;

end