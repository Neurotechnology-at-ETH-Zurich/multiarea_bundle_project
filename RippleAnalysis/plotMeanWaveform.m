function plotMeanWaveform(rippleEpisodes, indices)
% plotMeanWaveform plot the mean episode waveform and some random samples
%   rippleEpisodes  (n_ripples x n_samples)     array of ripples
%   indices         (:,)                        a vector of indices indicating which samples to plot

arguments
    rippleEpisodes
    indices (:,1) {mustBeInteger}
end

num_bg = 50; % the number of background samples to show
figure();
subset = rippleEpisodes(indices,:);
mean_waveform = mean(subset,1);
% check if we should show all episodes or collect a random subset
if numel(indices) <= num_bg
    plot(subset','Color',[0.5 0.5 0.5]);
else
    samples = randsample(indices,num_bg);
    subsubset = rippleEpisodes(samples,:);
    plot(subsubset','Color',[0.5 0.5 0.5]);
end
hold on;
plot(mean_waveform,'-b');
xlabel('frames');
ylabel('Potential');
title("Mean Waveform")
hold off;
end