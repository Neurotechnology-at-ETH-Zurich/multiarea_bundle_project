


 [ripples] = bz_FindRipples(lfp.data(2,:)',lfp.timestamps, 'durations',[20 100],'minDuration',5,'frequency',lfp.samplerate,'plotType',1,'passband',[180 200]); %Buzsaki uses 20ms for min duration
 [ripples_low] = bz_FindRipples(lfp.data(2,:)',lfp.timestamps, 'durations',[20 100],'minDuration',5,'frequency',lfp.samplerate,'plotType',1,'passband',[120 160]);

lfp.ripples=ripples;
lfp.ripples_low=ripples_low;

d=designfilt('bandpassfir','FilterOrder',600,'StopbandFrequency1',125,'PassbandFrequency1',130,'PassbandFrequency2',245,'StopbandFrequency2',250,'SampleRate',2000);
lfp.filtered = filtfilt(d,lfp.data(1,:));


index_segment_ripple_index=[];
separated_ripple_index=[];
for num_ripple=1:length (lfp.ripples_low.peaks(:,1))
    index_segment_ripple=intersect(find( lfp.ripples.timestamps(:,1) <=lfp.ripples_low.timestamps(num_ripple,1)),find(lfp.ripples_low.timestamps(num_ripple,1) <= lfp.ripples.timestamps(:,2)));
    index_segment_ripple_index=[index_segment_ripple_index; index_segment_ripple];
    if isempty(index_segment_ripple)
        separated_ripple_index=[separated_ripple_index; num_ripple];
    end
end

lfp.ripples_low_separated.peaks(:,1)=lfp.ripples_low.peaks(separated_ripple_index,1);
lfp.ripples_low_separated.timestamps(:,1:2)= lfp.ripples_low.timestamps(separated_ripple_index,1:2);

figure(3)
subplot(6,1,1)

plot(lfp.timestamps,lfp.data(1,:))
hold on; xline(ripples.timestamps(:,1),'r')
hold on; xline(ripples.timestamps(:,2),'b')

hold on; xline(lfp.ripples_low_separated.timestamps(:,1),'m')
hold on; xline(lfp.ripples_low_separated.timestamps(:,2),'k')
xlim([3589 3594])%xlim([2954 2957])
title('row data from iHP')
subplot(6,1,2)

plot(lfp.timestamps,lfp.filtered)
hold on; xline(ripples.timestamps(:,1),'r')
hold on; xline(ripples.timestamps(:,2),'b')

hold on; xline(lfp.ripples_low_separated.timestamps(:,1),'m')
hold on; xline(lfp.ripples_low_separated.timestamps(:,2),'k')
xlim([3589 3594])
title('filtered data from iHP')


%%plot single Ripple for the four structure:
 honnan=34 %(position.sample(Index_start,1)./position.amplifier_sample_rate);
%     [~,from]=min(abs((lfp.timestamps.*1000)-(3589.*1000)));
%      [~,to]=min(abs((lfp.timestamps.*1000)-(3594.*1000)));

from=find(ceil((lfp.timestamps).*1000)==(3589.*1000));
to=find(ceil((lfp.timestamps).*1000)==(3594.*1000));


    subplot(6,1,3,"replace")
    [cfs,f] = cwt(lfp.data(1,from(1):to(1)),'amor',1/DT,'VoicesPerOctave',32);
    helperCWTTimeFreqPlot(cfs,lfp.timestamps(from(1):to(1)),f,'surf','Intermediate hippocampus',[],'Frequency [Hz]');
    colormap('hot');
    caxis([0 20000]);
    ylim([0 250]);
    xlim([lfp.timestamps(from(1)) lfp.timestamps(to(1))])
%     hold on;
%     plot3([honnan honnan],[0 250],[200000 200000],'LineWidth',2,'Color','green','LineStyle','-');

       subplot(6,1,4,"replace")
    [cfs,f] = cwt(lfp.data(2,from(1):to(1)),'amor',1/DT,'VoicesPerOctave',32);
    helperCWTTimeFreqPlot(cfs,lfp.timestamps(from(1):to(1)),f,'surf','Dorsal hippocampus',[],'Frequency [Hz]');
    colormap('hot');
    caxis([0 20000]);
    ylim([0 250]);
    xlim([lfp.timestamps(from(1)) lfp.timestamps(to(1))])
 
%     hold on;
%     plot3([honnan honnan],[0 250],[200000 200000],'LineWidth',2,'Color','green','LineStyle','-');
%  
      subplot(6,1,5,"replace")
    [cfs,f] = cwt(lfp.data(3,from(1):to(1)),'amor',1/DT,'VoicesPerOctave',32);
    helperCWTTimeFreqPlot(cfs,lfp.timestamps(from(1):to(1)),f,'surf','Retrosplenial cortex ',[],'Frequency [Hz]');
    colormap('hot');
    caxis([0 20000]);
    ylim([0 250]);
    xlim([lfp.timestamps(from(1)) lfp.timestamps(to(1))])
%     hold on;
%     plot3([honnan honnan],[0 250],[200000 200000],'LineWidth',2,'Color','green','LineStyle','-');
%   
      subplot(6,1,6,"replace")
    [cfs,f] = cwt(lfp.data(4,from(1):to(1)),'amor',1/DT,'VoicesPerOctave',32);
    helperCWTTimeFreqPlot(cfs,lfp.timestamps(from(1):to(1)),f,'surf','Medial prefrontal cortex',[],'Frequency [Hz]');
    colormap('hot');
    caxis([0 20000]);
    ylim([0 250]);
    xlim([lfp.timestamps(from(1)) lfp.timestamps(to(1))])
%     hold on;
%     plot3([honnan honnan],[0 250],[200000 200000],'LineWidth',2,'Color','green','LineStyle','-');
hold off