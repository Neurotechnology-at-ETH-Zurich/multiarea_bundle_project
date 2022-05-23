function [x_pos, y_pos,trace_mat] = spikes_sortBin (spikes,spike_file,clean_threshold)
    
    load([spikes.basename,'.animal.behavior.mat']);
    for i=1:length(position.RostralXY_mm)-1
    d_rost(i) = pdist([position.RostralXY_mm(i,1:2);position.RostralXY_mm(i+1,1:2)],'euclidean');
    end
     
    position.RostralXY_mm(find(position.Rostral_XY(:,1)<clean_threshold(1)),1)=NaN;
    position.RostralXY_mm(find(position.Rostral_XY(:,1)<clean_threshold(1)),2)=NaN;
    position.RostralXY_mm(find(position.Rostral_XY(:,2)<clean_threshold(1)),1)=NaN;
    position.RostralXY_mm(find(position.Rostral_XY(:,2)<clean_threshold(1)),2)=NaN;
    position.RostralXY_mm([find(d_rost>clean_threshold(2));find(d_rost>clean_threshold(2))+1],1:2)=NaN;
     
     
    for i=1:length(position.CaudalXY_mm)-1
    d_caud(i) = pdist([position.CaudalXY_mm(i,1:2);position.CaudalXY_mm(i+1,1:2)],'euclidean');
    end
     
    position.CaudalXY_mm(find(position.CaudalXY_mm(:,1)<clean_threshold(1)),1)=NaN;
    position.CaudalXY_mm(find(position.CaudalXY_mm(:,1)<clean_threshold(1)),2)=NaN;
    position.CaudalXY_mm(find(position.CaudalXY_mm(:,2)<clean_threshold(1)),1)=NaN;
    position.CaudalXY_mm(find(position.CaudalXY_mm(:,2)<clean_threshold(1)),2)=NaN;
    position.CaudalXY_mm([find(d_caud>clean_threshold(3));find(d_caud>clean_threshold(3))+1],1:2)=NaN;
     
    %%filter the sample time by NaN location examples
    behavior = table2array(Position_Table);
    behavior(:,8:9) = position.RostralXY_mm;
    behavior(:,10:11) = position.CaudalXY_mm;
    
 
    % […] position.sample(find(~isnan(position.CaudalXY_mm(:,1))),1)./sample_rate;
    %only calculate firing rate batween valid timestamp. dT is the time between
    %the VALID time points....
    num_cell = spikes.numcells;
    num_sample = length(behavior(:,1));
    num_features = length(behavior(1,:));   
    
    %deleting rows with NaN
    
    ros_x = behavior(:,8);
    ros_y = behavior(:,9);
    
    
    

    x_pos = cell(1,num_cell);
    y_pos = cell(1,num_cell);
    
    ros_x(find(ros_x<0)) = 0;
    behavior(:,8) = ros_x;

    ros_y(find(ros_y<0)) = 0;
    behavior(:,9) = ros_y;
    
    
    behavior(:,[12,13]) = behavior(:,[12,13])/spikes.sr;
    trace_mat = horzcat(behavior(:,8),behavior(:,9),behavior(:,12));
    save([spikes.basename,'.trace.mat'],'trace_mat');
    disp('debug here');
    behavior_nan = find(isnan(behavior(:,8)));
    behavior(behavior_nan,:)=[];
    num_sample = length(behavior(:,1));
    
    %% bin
    edges = behavior(:,12);

    spikeTs = double(spike_file.spikeTimes)/spikes.sr;
    
    behavior = vertcat(behavior,zeros(1,num_features));
    
    
    %binning spikeTimes
    
    spikeTs_bin = discretize(spikeTs,edges);

    %getting rid of bad spikes by pointing them to the 0 values
    spikeTs_bin(isnan(spikeTs_bin)) = num_sample+1;

    spikes_xpos = behavior(:,8);
    spikes_xpos = spikes_xpos(spikeTs_bin);
    
    spikes_ypos = behavior(:,9);
    spikes_ypos = spikes_ypos(spikeTs_bin);
    

    max_spikes = length(spikeTs_bin);


    for i = 1:num_cell

        x_pos{i} = spikes_xpos(spikes.ids{i});
        x_pos{i} = x_pos{i}(find(x_pos{i} ~= 0));

        y_pos{i} = spikes_ypos(spikes.ids{i});
        y_pos{i} = y_pos{i}(find(y_pos{i} ~= 0));

    end
end