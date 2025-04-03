spikes_to_del = find(spikeClusters == -1);

fileID = fopen('combined_raw.jrc','r');
spikes_raw = fread(fileID, [41 11*nSpikes], 'int16');
spikes_raw = reshape(spikes_raw, [41 11 nSpikes]);
spikes_raw(:,:,spikes_to_del) = [];
fileID = fopen('combined_raw_clean.jrc','w');
fwrite(fileID, spikes_raw, 'int16');
clear spikes_raw

fileID = fopen('combined_filt.jrc','r');
spikes_filt = fread(fileID, [41 11*nSpikes], 'int16');
spikes_filt = reshape(spikes_filt, [41 11 nSpikes]);
spikes_filt(:,:,spikes_to_del) = [];
fileID = fopen('combined_filt_clean.jrc','w');
fwrite(fileID, spikes_filt, 'int16');
clear spikes_filt

fileID = fopen('combined_features.jrc','r');
features = fread(fileID, [2 21*nSpikes], 'single');
features = reshape(features, [3 14 nSpikes]);
features(:,:,spikes_to_del) = [];
fileID = fopen('combined_features_clean.jrc','w');
fwrite(fileID, features, 'single');
clear features

centerSites(spikes_to_del,:) = [];
initialClustering(spikes_to_del) = [];
ordRho(spikes_to_del) = [];
spikeSites(spikes_to_del) = [];
spikeSites2(spikes_to_del) = [];
spikeTimes(spikes_to_del) = [];
spikeAmps(spikes_to_del) = [];
spikeClusters(spikes_to_del) = [];
spikeDelta(spikes_to_del) = [];
spikeNeigh(spikes_to_del) = [];
spikePositions(spikes_to_del,:) = [];
spikeRho(spikes_to_del) = [];

spike_idx = 1:nSpikes;
spike_idx(spikes_to_del) = 0;
spike_idx = spike_idx';
spike_idx(spike_idx == 0) = [];
clusterCenters_new = zeros(size(clusterCenters));
for i=1:size(clusterCenters,1)
clusterCenters_new(i) = find(spike_idx == clusterCenters(i));
end
clusterCenters = clusterCenters_new;
clear clusterCenters_new

nClusters = size(clusterCenters,1);
for n=1:nClusters
   spikesByCluster_n = spikesByCluster{n};
   spikesByCluster_transposed = find(ismember(spike_idx, spikesByCluster_n));
   spikesByCluster{n} = spikesByCluster_transposed;
end
clear spikesByCluster_tranposed

nSites = size(spikesBySite,2);
for n=1:nSites
   spikesBySite_n = spikesBySite{n};
   spikesBySite2_n = spikesBySite2{n};
   spikesBySite_transposed = find(ismember(spike_idx, spikesBySite_n));
   spikesBySite2_transposed = find(ismember(spike_idx, spikesBySite2_n));
   spikesBySite{n} = spikesBySite_transposed;
   spikesBySite2{n} = spikesBySite2_transposed;
end
clear spikesBySite_transposed
clear spikesBySite2_transposed
clear nSites nClusters spike_idx n i 

nSpikes = nSpikes - size(spikes_to_del, 1);
featuresShape = [21,2,nSpikes];
filtShape = [41,11,nSpikes];
rawShape = [41,11,nSpikes];

clear spikes_to_del
save('combined_res_clean.mat')