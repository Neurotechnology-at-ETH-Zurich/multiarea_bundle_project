%% Load data from JRC file

load('amplifier_res.mat', 'meanWfLocal')
load('amplifier_res.mat', 'spikeClusters')
load('amplifier_res.mat', 'spikeTimes')
load('amplifier_res.mat', 'centerSites')

lfp.samplerate=2000;
bin=0.01; %parameter for PSTH (bar)
%% LOAD SINGLE Ch LFP 4 channels from the structres
extract_lfp=[12,81,220,205] %iHP, dHP, RSC, MPF (ch-1) in Neuroscope
for n_ch=1:length(extract_lfp)
    [lfp.data(n_ch,:)]=detrend(double(ImporterDAT('amplifier.lfp',256,extract_lfp(n_ch))))';
end
lfp.timestamps=(linspace(0,length(lfp.data)./lfp.samplerate,length(lfp.data)));
DT=1/lfp.samplerate;

%% selected cluID from SpikeCluster
cluID=25;

%2014 sec to 2016

%% upsample data

sample_for_waveform=squeeze(meanWfLocal(:,1,cluID))
upsampling_frequency = 44100;
evtWindowRaw = [-0.5, 1.5];

time_during_spike = (linspace(evtWindowRaw(1),evtWindowRaw(2),length(sample_for_waveform)))';
time_during_spike_upsampled=(linspace(evtWindowRaw(1),evtWindowRaw(2),(upsampling_frequency.*diff(evtWindowRaw ))./1000))';
%
upsampled_spike_wavform = interp1(time_during_spike,sample_for_waveform,time_during_spike_upsampled,'spline');

%plot(time_during_spike_upsampled,upsampled_spike_wavform)

res=spikeTimes(find(spikeClusters==25));

% filename = 'temp_2sec.wav';
% audiowrite(filename,upsampled_spike_wavform ,upsampling_frequency);
% wavewrite(upsampled_spike_wavform,44100,'test.WAV')

%% 2014
rajzolo=1; 
start_sec=2005 %secund -> segment of the video
stop_sec=2070

[~,Index_start]=min(abs(position.sample(:,1)-start_sec.*position.amplifier_sample_rate)); %index for the segment for position
[~,Index_stop]=min(abs(position.sample(:,1)-stop_sec.*position.amplifier_sample_rate));
row_data_wave=zeros(1,position.amplifier_sample_rate*(stop_sec-start_sec));
samples_of_spikes=(res(find(position.sample(Index_start) <= res & res <=position.sample(Index_stop))))-position.sample(Index_start); %samples of spike within a segment of the video

[~,hely]=min(sample_for_waveform);

for num_spikes_segment=1:length(samples_of_spikes)
    row_data_wave(samples_of_spikes(num_spikes_segment)-(hely-1):(samples_of_spikes(num_spikes_segment)-(hely-1)+length(sample_for_waveform)-1))=row_data_wave(samples_of_spikes(num_spikes_segment)-(hely-1):(samples_of_spikes(num_spikes_segment)-(hely-1)+length(sample_for_waveform)-1))+sample_for_waveform';
end

time_during_spike=(linspace(0,(stop_sec-start_sec),position.amplifier_sample_rate.*(stop_sec-start_sec)));
time_during_spike_upsampled=(linspace(0,(stop_sec-start_sec),(upsampling_frequency.*(stop_sec-start_sec))));
upsampled_spike_wavform = interp1(time_during_spike,row_data_wave,time_during_spike_upsampled,'spline');

filename = 'temp_1min.wav';
audiowrite(filename,upsampled_spike_wavform ,upsampling_frequency);
%wavewrite(upsampled_spike_wavform,44100,'test.WAV')
% [AUDIO,Fs] = audioread('temp_2sec.wav');

%% spikes/bin for magnitude


for index_frames=1:length(position.sample)-1
    spikes_per_bin=res(find(position.sample(index_frames,1)<= res & res < position.sample(index_frames+1,1)));
    position.angle(index_frames,3)=numel(spikes_per_bin);
end

segment_angles_spikes=[(Index_start:Index_stop)', position.angle(Index_start:Index_stop,:)]; %original index, angle, rad, spikes/frame (relative index is the row number)
cm = jet(max(unique(segment_angles_spikes(:,4))).*10);   %color for arrows set by the maximum firing/frame *10

%temp_vector_for_polar=segment_angles_spikes(find(segment_angles_spikes(:,4)>0),1:4);
%total_angle=1;
% for ploar_vect_length=1:length(temp_vector_for_polar)
%     multiplicator=temp_vector_for_polar(ploar_vect_length,4)
%     for accumulate_angle=1:multiplicator
% pol_input(total_angle)=temp_vector_for_polar(ploar_vect_length,3);
% total_angle=total_angle+1;
%     end
% end
%
% polarhistogram(pol_input+(pi/2),100,'FaceColor','white','FaceAlpha',.8)
% set(gca,'Color','k');
% set(gcf,'color','k');
% ax = gca % Get handle to current axes.
% ax.GridAlpha = 0.9;  % Make grid lines less transparent.
% ax.GridColor = [0.1, 0.7, 0.2];
% ax.Color = 'k';
% ax.RColor =  [0.1, 0.7, 0.2];
% ax.FontSize = 14
% ax. ThetaColor = [0.1, 0.7, 0.2];
% ax.RLim =[0 140]

% honnan=floor(position.sample(Index_start,1)./position.amplifier_sample_rate);
% from=find(ceil((lfp.timestamps).*1000)==(honnan-0.5).*1000);
% to=find(ceil((lfp.timestamps).*1000)==(honnan+0.5).*1000);

% [cfs,f] = cwt(lfp.data(1,from(1):to(1)),'amor',1/DT,'VoicesPerOctave',32);
% helperCWTTimeFreqPlot(cfs,lfp.timestamps(from(1):to(1)),f,'surf','LFP (CWT)','seconds','Hz');
% ylim([0 250]);
% colormap('hot');
% caxis([0 8000]);
% hold on;
% plot3([honnan honnan],[0 250],[200000 200000],'LineWidth',4,'Color','green','LineStyle','--');

clear color_accumulation    nyil_accumulation nyil_accumulation_magnitude    pol_input
pol_input=[]
v = VideoWriter(['test3.avi']);
v.FrameRate = position.video_sample_rate;
open(v);
segment_counter=1;
total_angle=1
while Index_start<= Index_stop

    honnan=(position.sample(Index_start,1)./position.amplifier_sample_rate);
    [~,from]=min(abs((lfp.timestamps.*1000)-((honnan-1).*1000)));
     [~,to]=min(abs((lfp.timestamps.*1000)-((honnan+1).*1000)));

    subplot(5,2,[1],"replace")
    [cfs,f] = cwt(lfp.data(1,from(1):to(1)),'amor',1/DT,'VoicesPerOctave',32);
    helperCWTTimeFreqPlot(cfs,lfp.timestamps(from(1):to(1)),f,'surf','Intermediate hippocampus',[],'Frequency [Hz]');
    colormap('hot');
    caxis([0 8000]);
    ylim([0 200]);
    xlim([lfp.timestamps(from(1)) lfp.timestamps(from(1)+4001)])
%     hold on;
%     plot3([honnan honnan],[0 250],[200000 200000],'LineWidth',2,'Color','green','LineStyle','-');

    subplot(5,2,[3],"replace")
    [cfs,f] = cwt(lfp.data(2,from(1):to(1)),'amor',1/DT,'VoicesPerOctave',32);
    helperCWTTimeFreqPlot(cfs,lfp.timestamps(from(1):to(1)),f,'surf','Dorsal hippocampus',[],'Frequency [Hz]');
    colormap('hot');
    caxis([0 8000]);
    ylim([0 200]);
    xlim([lfp.timestamps(from(1)) lfp.timestamps(from(1)+4001)])
 
%     hold on;
%     plot3([honnan honnan],[0 250],[200000 200000],'LineWidth',2,'Color','green','LineStyle','-');
%  
    subplot(5,2,[5],"replace")
    [cfs,f] = cwt(lfp.data(3,from(1):to(1)),'amor',1/DT,'VoicesPerOctave',32);
    helperCWTTimeFreqPlot(cfs,lfp.timestamps(from(1):to(1)),f,'surf','Retrosplenial cortex ',[],'Frequency [Hz]');
    colormap('hot');
    caxis([0 8000]);
    ylim([0 200]);
    xlim([lfp.timestamps(from(1)) lfp.timestamps(from(1)+4001)])
%     hold on;
%     plot3([honnan honnan],[0 250],[200000 200000],'LineWidth',2,'Color','green','LineStyle','-');
%   
    subplot(5,2,[7],"replace")
    [cfs,f] = cwt(lfp.data(4,from(1):to(1)),'amor',1/DT,'VoicesPerOctave',32);
    helperCWTTimeFreqPlot(cfs,lfp.timestamps(from(1):to(1)),f,'surf','Medial prefrontal cortex',[],'Frequency [Hz]');
    colormap('hot');
    caxis([0 8000]);
    ylim([0 200]);
    xlim([lfp.timestamps(from(1)) lfp.timestamps(from(1)+4001)])
%     hold on;
%     plot3([honnan honnan],[0 250],[200000 200000],'LineWidth',2,'Color','green','LineStyle','-');
hold off

subplot(5,2,[9])
if rajzolo==1
spikeTimes=double(spikeTimes);
spikeClusters=double(spikeClusters);

Neuron.Times=double(spikeTimes(spikeClusters>0))./20000;
Neuron.Clusters=double(spikeClusters(spikeClusters>0));

edges=[min(lfp.timestamps):bin:max(lfp.timestamps)]; 

for neurons=1:max(Neuron.Clusters)
    if clusterSites(neurons)<=64  && neurons ~=cluID
        col='c';
        marker='.';
    elseif clusterSites(neurons)>64 && clusterSites(neurons)<=128  && neurons ~=cluID
        col='r';
        marker='.';
    elseif clusterSites(neurons)>128 && clusterSites(neurons)<=192  && neurons ~=cluID
        col='g';
        marker='.';
    elseif clusterSites(neurons)>192 && neurons ~=cluID
        col='w';
        marker='.';
    elseif  neurons==cluID
        marker='.';
        col='y';
    end
   scatter(Neuron.Times(Neuron.Clusters==neurons),repmat(neurons,1,length(find(Neuron.Clusters==neurons))), marker,'MarkerEdgeColor',col,'MarkerFaceColor',col)
    bined_spikes(neurons,:)=histcounts(Neuron.Times(Neuron.Clusters==neurons),edges);
  

  % xlim([2970 2972])
hold on
end
xlim([lfp.timestamps(from(1)) lfp.timestamps(from(1)+4001)]);
set(gca, 'Color','k', 'XColor','w', 'YColor','w')
rajzolo=0;
else
    xlim([lfp.timestamps(from(1)) lfp.timestamps(from(1)+4001)]);
end

    obj.CurrentTime= (Index_start-1)./position.video_sample_rate
    refImage = readFrame(obj);

    % scatter(position.Rostral_XY(Index_start,1),position.Rostral_XY(Index_start,2),'.g');
    % scatter(position.Caudal_XY(Index_start,1),position.Caudal_XY(Index_start,2),'.r');

    Y1 =  position.Caudal_XY(Index_start,2);
    X1=  position.Caudal_XY(Index_start,1);
    Y2 = position.Rostral_XY(Index_start,2);
    X2 = position.Rostral_XY(Index_start,1);
    delta_X = X2 - X1;
    delta_Y = Y2 - Y1;

    % theta = atan2(Y1-Y2 ,X1-X2);
    % rad2deg(theta)

    %magnitude = sqrt(((delta_X).^2)+((delta_Y).^2));
    % create colormap
    if  0<segment_angles_spikes(segment_counter,4)
        ind = segment_angles_spikes(segment_counter,4)*10;
    else
        ind = 1;    % convert magnitude to index
    end
    %quiver(X1,Y1,delta_X,delta_Y,'color',cm(ind,:),'LineWidth',1.2)
        subplot(5,2,[6 8 10])
        set(gcf,'color','k');
    if abs(delta_X)<150 && abs(delta_Y)<150
        nyil_accumulation(segment_counter,:)=[X1,Y1, 0];
        nyil_accumulation_magnitude(segment_counter,:)=[delta_X,delta_Y,0];
        color_accumulation(segment_counter,:)=cm(ind,:);

        temp_vector_for_polar(segment_counter,:)=segment_angles_spikes(segment_counter,:);
        multiplicator=temp_vector_for_polar(segment_counter,4);
        for accumulate_angle=1:multiplicator
            pol_input(total_angle)=temp_vector_for_polar(segment_counter,3);
            total_angle=total_angle+1;
        end
        if segment_counter<=10
            hold off
            imshow(refImage)
            hold on;
            qHandle = quiver3D(nyil_accumulation, nyil_accumulation_magnitude, color_accumulation);
            lighting phong;
            camlight head;
            segment_counter=segment_counter+1
        elseif  segment_counter>10
            hold off
            imshow(refImage)
            hold on;
            qHandle = quiver3D(nyil_accumulation(segment_counter-10:segment_counter,:), nyil_accumulation_magnitude(segment_counter-10:segment_counter,:), color_accumulation(segment_counter-10:segment_counter,:));
            lighting phong;
            camlight head;
            segment_counter=segment_counter+1
        end

        subplot(5,2,[2 4],"replace")
        polarhistogram(pol_input-(pi/2),50,'FaceColor','white','FaceAlpha',1)
        %  polarhistogram(pol_input-((3*pi)/2),50,'FaceColor','white','FaceAlpha',1)
        ax = gca % Get handle to current axes.
        ax.GridAlpha = 0.9;  % Make grid lines less transparent.
        ax.GridColor = [0.1, 0.7, 0.2];
        ax.Color = 'k';
        ax.RColor =  [0.1, 0.7, 0.2];
        ax.FontSize = 14
        ax. ThetaColor = [0.1, 0.7, 0.2];
%        ax.RLim =[0 140]
    else
          imshow(refImage)
    end


   
    set(gcf,'PaperPositionMode','manual','PaperUnits','points','PaperSize',[1920 1200],'Position',[1,1,1920,1200]);
    videoFrame = getframe(gcf);
    writeVideo(v, videoFrame);
    Index_start=Index_start+1;
end
close(v)
%
% %%
% frames=Index_stop-Index_start
% numAudio = size(upsampled_spike_wavform,2);
% numRep = floor(numAudio/(frames));
% numDiff = numAudio - numRep*(frames); % mismatch
%  upsampled_spike_wavform = upsampled_spike_wavform';
% if numDiff~=0
%     % if length(frames) does not evenly divide nAudioSamples, then
%     % subsample audio to match numRep*length(frames)
%     selector = round(linspace(1, numAudio, numRep*(frames)));
%     subSignal = upsampled_spike_wavform (selector, :);
% end
% assert(numRep*(frames) == size(subSignal,1));
%
% shotPath=pwd
% videoFWriter = vision.VideoFileWriter(fullfile(shotPath, 'avclip','av_clip.avi'), ...
%                                       'AudioInputPort', true);
% for i = 1:length(frames)
%    fprintf('Frame: %d/%d\n', i, length(frames));
%    step(videoFWriter, frames{i}, subSignal(numRep*(i-1)+1:numRep*i,:));
% end
%
% while Index_start<=Index_stop
% obj.CurrentTime= (Index_start-1)./position.video_sample_rate
% refImage = readFrame(obj);
% imshow(refImage)
% hold on;
% scatter(position.Rostral_XY(Index_start,1),position.Rostral_XY(Index_start,2),'*g');
% scatter(position.Caudal_XY(Index_start,1),position.Caudal_XY(Index_start,2),'*r');
% set(gcf,'PaperPositionMode','manual','PaperUnits','points','PaperSize',[1920 1200],'Position',[1,1,1920,1200]);
% videoFrame = getframe(gcf);
% writeVideo(v, frame);
% Index_start=Index_start+1;
% end
%
% release(videoFWriter);
%
%
%
%
% %       alpha(0.3)
% frame = getframe(gcf);
% writeVideo(v, frame);
%
%
% videoFReader = vision.VideoFileReader('video1.avi');
% [AUDIO,Fs] = audioread('audio1.wav');
% videoFWriter = VideoFileWriter('newvideo.avi','AudioInputPort',true);
% videoFrame = step(videoFReader);
% step(videoFWriter,videoFrame,AUDIO);
% release(videoFReader);
% release(videoFWriter);