function alignedCenterFrames  = alignRipples(signal,rippleCenterSeeds,eventWindowFrames)
% alignRipples perform cross correlation based alignment of ripples
% signal (time,) signal containing the events of interest
% rippleCenterSeeds (ripples, ) the inital center frames
% eventWindow (1,) the duration of the windows around the center frames
% used for alignment
%
% alignedCenterFrames (ripples,) the center frames of the aligned
% eventWindows

% start by using the provided seeds as central frames
rippleEpisodes.centralFrames = rippleCenterSeeds;

% extract centered windows around the seeds for each ripple
eventWindowHalfFrames = round(eventWindowFrames/2);
num_ripples = size(rippleCenterSeeds,1);
for i = 1:num_ripples
    startFrame = rippleEpisodes.centralFrames(i) - eventWindowHalfFrames;
    endFrame = rippleEpisodes.centralFrames(i) + eventWindowHalfFrames;
    rippleEpisodes.episodes(i,:) = signal(startFrame:endFrame);
end

% Iteratively refine the center frames by maximizing cross correlation to a
% signal template constructed by averaging all ripples that are aligned
% restrict shift to +/- half a ripple cycle
% ripples ~= 200 Hz  = 10 frame cycle at 2000Hz sampling rate
maxlag = 2; % < half cycle of ripple signal in frames for stability


f = waitbar(0, 'Aligning ripples...');
for iteration = 1:20

    % calculate alignment template (time x channels)
    mean_waveform = mean(rippleEpisodes.episodes,1);

    % emphasize alignment in the central region of the ripple by tampering the
    % signal using a hamming window
    template = mean_waveform(:) .* hamming(numel(mean_waveform));

    % greedily align signals to template by shifting them up to maxlag
    for i = 1:num_ripples
        % find signal offset by maximizing cross correlation
        r = xcorr(template,rippleEpisodes.episodes(i,:),maxlag);
        [maxr,d] = max(r); % get highest cross correlation and its associated delay
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
    

    % update progress bar
    waitbar(i/num_ripples, f, sprintf('Progress: %d %%', floor(i/num_ripples*100)));
end
close(f)

% return the aligned center frames
alignedCenterFrames = rippleEpisodes.centralFrames;

end