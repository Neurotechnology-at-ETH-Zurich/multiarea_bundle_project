%% Path for the Sessions from TBY_37 rat for Assmebly detection
%CONSTANT
path_counter=1;
%% kernel for smoothing
sigma = .005 ;                            %Standard deviation of the kernel 5-15 ms
edges_norm=[-0.5*sigma:.001:0.5*sigma];     %Time ranges form -X0*st. dev. to X1*st. dev.
kernel = normpdf(edges_norm,0,sigma);
kernel = kernel*.001;
edges=[-0.5:0.005:0.5];
edges_CCA=[-0.5:0.025:0.5];
%edges_CCA_ONE_MS=[-0.5:0.01:0.5];
sample_rate=20000;

sigma_for_trial=0.005;
bin_trial=0.001;
edges_trial=[-0.5:0.015:0.5];
edges_norm_for_trial=[-0.5*sigma_for_trial:bin_trial:0.5*sigma_for_trial];
kernel_trial = normpdf(edges_norm_for_trial,0,sigma_for_trial);
kernel_trial=kernel_trial*bin_trial;

%% CONSTANTS
sample_rate=20000;
lfp.samplerate=2000;
bin=0.01; %parameter for PSTH (bar)
CellID_counter=1;
%   edges_CCA=[-0.5:0.025:0.5];
%% LOAD LFP and Spikes

%you have to be in the folder where selected session are.
base_folder=pwd;
listing = dir('F:\Dropbox (Personal)\ETH_DATA\rTBY34');

% getting the folder with the data, eliminate noise from path
for i=1:numel(listing)
    k = strfind(listing(i).name,'freely_behav');
    if ~isempty(k)
        path(path_counter).folders=fullfile(base_folder,listing(i).name);
        l=strfind(listing(i).name,'_');
        path(path_counter).sessionID=str2num(listing(i).name(1:l(1)-1));
        path_counter=path_counter+1;
    end

end

% ordered of the session by session ID
[~,sorting_session]=sort([path.sessionID]);
path=path(sorting_session);

%clear up memory
clear k l i listing path_counter sorting_session

%% Grabbing the files: JRC, Ripples events previously corated, features table, lfp etc

for nume_files=1:numel(path)

    Path_temp=char(path(nume_files).folders);
    %     disp(['Load LFP from: '  char(fullfile(Path_temp,'amplifier.lfp'))]);
    %
    %     lfp.data(nume_files)={detrend(double(ImporterDAT(char(fullfile(Path_temp,'amplifier.lfp')),256,selected_channel_fromHP(nume_files))))'};                     % Row LFPs filtered 0.1-1000Hz during saving [int16]
    %     lfp.timestamps(nume_files)={(linspace(0,length([lfp.data{nume_files}(:,:)])./lfp.samplerate,length([lfp.data{nume_files}(:,:)])))};   %time stamps [double]

    listing =dir(Path_temp);

    %clearin results of dir, only folders

    % Get ripples
    for l=1:numel(listing)
        if ~isempty(strfind(listing(l).name,'auto_features'))
            disp(['Load RIPPLES data from: '  char(fullfile(Path_temp,(listing(l).name)))]);
            Ripples(nume_files)=load(fullfile(Path_temp,(listing(l).name)),'featureTable','lfp','lfp_main_channel_dHP','lfp_main_channel_iHP','ripples');
        end
    end

    %clearin results of dir, only folders
    dirFlags =[listing.isdir];

    listing=listing(dirFlags);
    listing=listing(~ismember({listing.name},{'.','..'}));

    % Features_Ripples=
    % dir ['Path_temp\*features*.m]

    for i=1:numel(listing)
        k = strfind(listing(i).name,'RSC');
        if ~isempty(k) & isempty(strfind(listing(i).name,'except'))
            disp(['Load Spikes from: '  char(fullfile( Path_temp,listing(i).name,'amplifier_res.mat'))]);
            Spike_RSC.Sites(nume_files)=load(char(fullfile( Path_temp,listing(i).name,'amplifier_res.mat')), 'clusterSites');
            Spike_RSC.Cluster(nume_files)=load(char(fullfile( Path_temp,listing(i).name,'amplifier_res.mat')), 'spikeClusters');
            Spike_RSC.Times(nume_files)=load(char(fullfile( Path_temp,listing(i).name,'amplifier_res.mat')), 'spikeTimes');
            Spike_RSC.clusterNotes(nume_files)=load(char(fullfile( Path_temp,listing(i).name,'amplifier_res.mat')), 'clusterNotes');
            Spike_RSC.meanWfLocalRaw(nume_files)=load(char(fullfile( Path_temp,listing(i).name,'amplifier_res.mat')), 'meanWfLocalRaw');
        end
        if ~isempty(k) & ~isempty(strfind(listing(i).name,'except'))
            disp(['Load Spikes from: '  char(fullfile( Path_temp,listing(i).name,'amplifier_res.mat'))]);
            Spike_Excluded.Sites(nume_files)=load(char(fullfile( Path_temp,listing(i).name,'amplifier_res.mat')), 'clusterSites');
            Spike_Excluded.Cluster(nume_files)=load(char(fullfile( Path_temp,listing(i).name,'amplifier_res.mat')), 'spikeClusters');
            Spike_Excluded.Times(nume_files)=load(char(fullfile( Path_temp,listing(i).name,'amplifier_res.mat')), 'spikeTimes');
            Spike_Excluded.clusterNotes(nume_files)=load(char(fullfile( Path_temp,listing(i).name,'amplifier_res.mat')), 'clusterNotes');
            Spike_Excluded.meanWfLocalRaw(nume_files)=load(char(fullfile( Path_temp,listing(i).name,'amplifier_res.mat')), 'meanWfLocalRaw');
        end
    end

    % temp_hely=find((Spike_Excluded.Cluster.spikeClusters==45));
    % unique(Spike_Excluded.Times(nume_files).spikeTimes(temp_hely))
    %
    % unique(Spike_Excluded.Sites(nume_files).spikeSites(temp_hely))
    % merging files RSC Excluded RSC
    Spike.Sites(nume_files).spikeSites=[Spike_Excluded.Sites(nume_files).clusterSites Spike_RSC.Sites(nume_files).clusterSites];
    Spike.Times(nume_files).spikeTimes=double([Spike_Excluded.Times(nume_files).spikeTimes; Spike_RSC.Times(nume_files).spikeTimes]);
    Spike_RSC.Cluster(nume_files).spikeClusters(find(Spike_RSC.Cluster(nume_files).spikeClusters>0))= Spike_RSC.Cluster(nume_files).spikeClusters(find(Spike_RSC.Cluster(nume_files).spikeClusters>0))+max(unique(Spike_Excluded.Cluster(nume_files).spikeClusters));
    Spike.Cluster(nume_files).spikeCluster=[Spike_Excluded.Cluster(nume_files).spikeClusters; Spike_RSC.Cluster(nume_files).spikeClusters];
    Spike.clusterNotes(nume_files).clusterNotes=[Spike_Excluded.clusterNotes(nume_files).clusterNotes; Spike_RSC.clusterNotes(nume_files).clusterNotes];
    Spike.meanWfLocalRaw(nume_files).meanWfLocalRaw=[squeeze(Spike_Excluded.meanWfLocalRaw(nume_files).meanWfLocalRaw(:,1,:)) squeeze(Spike_RSC.meanWfLocalRaw(nume_files).meanWfLocalRaw(:,1,:))];


    %fetching variables for the Cells matrix

    spikeSites=Spike.Sites(nume_files).spikeSites;
    spikeTimes=double(Spike.Times(nume_files).spikeTimes);
    spikeClusters=Spike.Cluster(nume_files).spikeCluster;
    waveform=Spike.meanWfLocalRaw(nume_files).meanWfLocalRaw;
    ripples.eventID=Ripples(nume_files).ripples.eventID;
    featureTable=Ripples(nume_files).featureTable;

    Recording_CellID=(unique(spikeClusters(spikeClusters>0)));
    PutitativeStructure={};

    for CellID=1:size(Recording_CellID,1)
        animal_data.CellID=Recording_CellID(CellID);

        %         if num_lfp==1 || num_lfp==2 || num_lfp==3
        %             if unique(unique(spikeSites(spikeClusters==CellID))<=64)
        %                 PutitativeStructure(CellID)={'mPFC'};
        %             elseif unique(unique(spikeSites(spikeClusters==CellID))>64 & unique(spikeSites(spikeClusters==CellID))<=128)
        %                 PutitativeStructure(CellID)={'RSC'};
        %             elseif unique(unique(spikeSites(spikeClusters==CellID))>128 & unique(spikeSites(spikeClusters==CellID))<=192)
        %                 PutitativeStructure( CellID)={'dHP'};
        %             elseif unique(unique(spikeSites(spikeClusters==CellID))>192)
        %                 PutitativeStructure( CellID)={'iHP'};
        %             end
        %         elseif num_lfp==2
        %             if unique(unique(spikeSites(spikeClusters==CellID))<=64)
        %                 PutitativeStructure(CellID)={'iHP'};
        %             elseif unique(unique(spikeSites(spikeClusters==CellID))>64 & unique(spikeSites(spikeClusters==CellID))<=128)
        %                 PutitativeStructure(CellID)={'dHP'};
        %             elseif unique(unique(spikeSites(spikeClusters==CellID))>128 & unique(spikeSites(spikeClusters==CellID))<=192)
        %                 PutitativeStructure( CellID)={'RSC'};
        %             elseif unique(unique(spikeSites(spikeClusters==CellID))>192)
        %                 PutitativeStructure( CellID)={'mPFC'};
        %             end
        %         elseif num_lfp==3
        if spikeSites(CellID)<=64
            PutitativeStructure(CellID)={'iHP'};
        elseif spikeSites(CellID)>64 & spikeSites(CellID)<=128
            PutitativeStructure(CellID)={'dHP'};
        elseif  spikeSites(CellID)>128 &  spikeSites(CellID)<=192
            PutitativeStructure( CellID)={'RSC'};
%         elseif  spikeSites(CellID)>192
%             PutitativeStructure( CellID)={'mPFC'};
                    elseif  spikeSites(CellID)>192  &  spikeSites(CellID)<=202
                        PutitativeStructure( CellID)={'IL'};
                    elseif   spikeSites(CellID)>202  &  spikeSites(CellID)<=220
                        PutitativeStructure( CellID)={'PrL'};
                    elseif   spikeSites(CellID)>220
                        PutitativeStructure( CellID)={'Cg1'};
        end
        %         end
        animal_data.PutitativeStructure=PutitativeStructure;


        %         time_on=(lfp.ripples{num_lfp}.timestamps(:,1))-0.5;
        %         time_off=( lfp.ripples{num_lfp}.timestamps(:,1))+0.5;%(ripples.timestamps(:,2));
        %         center= lfp.ripples{num_lfp}.timestamps(:,1); %(ripples.peaks); %center of peak

        time_on=featureTable.("center Timepoint")-0.5;
        time_off=featureTable.("center Timepoint")+0.5;%(ripples.timestamps(:,2));
        center= featureTable.("center Timepoint"); %(ripples.peaks); %center of peak


        center_rand=Ripples(nume_files).lfp.timestamps(randi(length(Ripples(nume_files).lfp.timestamps),1,length(time_on)))';
        while ~isempty(intersect(center_rand, center)) % Avoid intersection with the random time points (at center location!!!)
            center_rand=Ripples(nume_files).lfp.timestamps(randi(length(Ripples(nume_files).lfp.timestamps),1,length(time_on)));
        end

        time_on_rand=center_rand-0.5;
        time_off_rand=center_rand+0.5;


        fprintf(' in CellID %d |\n',CellID);

        %  declare variables, and counters
        trial_ripple_ON=1;
        SpikeSecund=(spikeTimes(find(spikeClusters==CellID)))./sample_rate;

        if length(SpikeSecund)>100
            spikes_trials=zeros(length(center),length(edges_CCA));
            spikes_trials_random=zeros(length( center),length(edges_CCA));
           % spikes_trials_one_ms=zeros(length( center),length(edges_CCA_ONE_MS));
            spikes_trials_assembly=zeros(length(center),length(edges_trial));
            spikes_trials_assembly_random=zeros(length(center),length(edges_trial));
            spikes_trials_assembly_Gauss=zeros(length(center),length(edges_trial));
            spikes_trials_assembly_Gauss_rand=zeros(length(center),length(edges_trial));
            for trial=1:length(center)
                fprintf_r('Ripples %i', trial)


                if ~isempty(SpikeSecund (find(time_on(trial) <= SpikeSecund   & SpikeSecund   <= time_off(trial))))
                    spike(trial_ripple_ON,:).time=(SpikeSecund(find(time_on(trial) <= SpikeSecund & SpikeSecund <=time_off(trial))))-(center(trial));
                    spike(trial_ripple_ON,:).counts_pre=[length(find(-0.5<= spike(trial_ripple_ON,:).time & spike(trial_ripple_ON,:).time<-0.3)  ) ];
                    spike(trial_ripple_ON,:).counts_post=[length(find(-0.1< spike(trial_ripple_ON,:).time & spike(trial_ripple_ON,:).time <=0.1))];
                    spike(trial_ripple_ON,:).counts_pre_test=[length(find(-0.5<= spike(trial_ripple_ON,:).time & spike(trial_ripple_ON,:).time<-0.3))];
                    spike(trial_ripple_ON,:).counts_post_test=[length(find(-0.1< spike(trial_ripple_ON,:).time & spike(trial_ripple_ON,:).time <=0.1))];
                    spike(trial_ripple_ON,:).reipple_ID=[trial, nume_files];
                    spike(trial_ripple_ON,:).ripple_classes=[char(ripples.eventID(trial))];
                    %                 spikes_gauss(trial_ripple_ON,:)=conv(histc(spike(trial_ripple_ON,:).time,edges)',kernel);
                    spikes_gauss(trial_ripple_ON,:)=(histc(spike(trial_ripple_ON,:).time,edges)');

                    spikes_trials(trial,:)=(histc(spike(trial_ripple_ON,:).time,edges_CCA)');
                   % spikes_trials_one_ms(trial,:)=(histc(spike(trial_ripple_ON,:).time,edges_CCA_ONE_MS)');

                    spikes_trials_assembly(trial,:)=(histc(spike(trial_ripple_ON,:).time,edges_trial)');
                    s=conv(spikes_trials_assembly(trial,:),kernel_trial);
                    center_gauss= ceil(length(edges_trial)/2);
                    spikes_trials_assembly_Gauss(trial,:)=s(ceil(length(s)/2)-( center_gauss-1):ceil(length(s)/2)+( center_gauss-1));
                    trial_ripple_ON=trial_ripple_ON+1;

                end

                if ~isempty(SpikeSecund (find(time_on_rand(trial) <= SpikeSecund   & SpikeSecund   <= time_off_rand(trial))))
                     spikes_trials_random(trial,:)=(histc((SpikeSecund(find(time_on_rand(trial) <= SpikeSecund & SpikeSecund <=time_off_rand(trial))))-(center_rand(trial)),edges_CCA))';
                    spikes_trials_assembly_random(trial,:)=(histc((SpikeSecund(find(time_on_rand(trial) <= SpikeSecund & SpikeSecund <=time_off_rand(trial))))-(center_rand(trial)),edges_trial))';
                    s_rand=conv(spikes_trials_assembly_random(trial,:),kernel_trial);
                    center_gauss= ceil(length(edges_trial)/2);
                    spikes_trials_assembly_Gauss_rand(trial,:)=s_rand(ceil(length(s_rand)/2)-(center_gauss-1):ceil(length(s_rand)/2)+(center_gauss-1));
          
                end


            end

            temp_conv_bzHigh_dHP=conv(zscore(mean(spikes_gauss(ismember({spike.ripple_classes},{'bzHigh_dHP'}),1:end-1),1)),kernel);
            temp_conv_bzHigh_iHP=conv(zscore(mean(spikes_gauss(ismember({spike.ripple_classes},{'bzHigh_iHP'}),1:end-1),1)),kernel);
            temp_conv_bzLow_dHP=conv(zscore(mean(spikes_gauss(ismember({spike.ripple_classes},{'bzLow_dHP'}),1:end-1),1)),kernel);
            temp_conv_bzLow_iHP=conv(zscore(mean(spikes_gauss(ismember({spike.ripple_classes},{'bzLow_iHP'}),1:end-1),1)),kernel);
            Cells(CellID_counter).spikes_kernel_mean=[temp_conv_bzHigh_dHP(3:end-3); temp_conv_bzHigh_iHP(3:end-3);   temp_conv_bzLow_dHP(3:end-3); temp_conv_bzLow_iHP(3:end-3)];
            Cells(CellID_counter).spikeTimesSec={spike.time};
            Cells(CellID_counter).ripple_classes=  {spike.ripple_classes}';
            temp=[spike.reipple_ID];
            temp_rippleID=[temp(1:2:end); temp(2:2:end)]';

            Cells(CellID_counter).Ripple_ID=  temp_rippleID;
            Cells(CellID_counter).individual_trials=spikes_trials;
            Cells(CellID_counter).individual_trials_random= spikes_trials_random;
       %     Cells(CellID_counter).individual_trials_one_ms= spikes_trials_one_ms;
            Cells(CellID_counter).individual_trials_count_assembly=  spikes_trials_assembly;
            Cells(CellID_counter).individual_trials_gauss_assembly=  spikes_trials_assembly_Gauss;
            Cells(CellID_counter).individual_trials_count_assembly_rand=    spikes_trials_assembly_random;
            Cells(CellID_counter).individual_trials_gauss_assembly_rand=   spikes_trials_assembly_Gauss_rand;

            [p,h,stats] = ranksum([spike.counts_pre_test],[spike.counts_post_test]);
            try
                Cells(CellID_counter).h_p_zval_ranksum=[h,p,stats.zval,stats.ranksum];
            catch
                Cells(CellID_counter).h_p_zval_ranksum=[h,p,NaN,stats.ranksum];
            end
            Cells(CellID_counter).modulation=(sum([spike.counts_post])-sum([spike.counts_pre]))./(sum([spike.counts_post])+sum([spike.counts_pre]));
            Cells(CellID_counter).Time=edges(1:end-1);
            Cells(CellID_counter).Structre=char(animal_data.PutitativeStructure(CellID));
            Cells(CellID_counter).CellID_original=CellID;
            Cells(CellID_counter).spikeSites=spikeSites(CellID);
            Cells(CellID_counter).clusterNotes=Spike.clusterNotes(nume_files).clusterNotes{CellID}(:,:);
            Cells(CellID_counter).wavform= waveform(:,CellID);
            Cells(CellID_counter).session=(path(nume_files).sessionID);
            Cells(CellID_counter).ratID='rTBY37';



            CellID_counter=CellID_counter+1;
        end
          clear spikes_gauss spike ic trial temp  temp_rippleID spikes_trials temp_conv_bzHigh_dHP temp_conv_bzHigh_iHP  temp_conv_bzLow_dHP temp_conv_bzLow_iHP  spikes_trials_one_ms
        fprintf_r('reset')

    end
    clear time_on time_off  center featureTable waveform  ripples time_on_rand time_off_rand center_rand

end

%% correcting clusterNotes where label is missing
emptyCells = find(cellfun(@isempty,{Cells.clusterNotes}));
for num_error=1:length(emptyCells)
    Cells(emptyCells(num_error)).clusterNotes='nothing';
end

[session_ID,~,session_counter]=unique([Cells.session]);

%% Assembly Detections

start_dtop_rip=[find(round(edges_CCA,2)==-0.13) find(round(edges_CCA,2)==0.13)];


len=length(edges_CCA)*2; %Random segments are concataneded to the actual segment to have the same assembly detection, therfore x2
opts.threshold.method = 'MarcenkoPastur';
opts.Patterns.method = 'ICA';
opts.Patterns.number_of_iterations = 10000;

disp('Single units from mPFC, RSC, dHP, iHP -> concataneted activity during ripples bin 25ms')

structures_unique=unique({Cells.Structre});

for nume_files=1:numel(path)


    for structures_NUM=1:length(structures_unique)
        darab_neuron=1;
        for num_cell=intersect(intersect(find(ismember({Cells.clusterNotes},{'single'})), find(ismember({Cells.Structre},structures_unique(structures_NUM)))),find(ismember([Cells.session],session_ID(nume_files))))
            disp(['cell in ' char(structures_unique{structures_NUM}) ' ID:' num2str(num_cell)])
            row_selected_ripples=ismember(Ripples(session_counter(num_cell)).ripples.eventID,{'bzHigh_dHP','bzLow_dHP','bzLow_iHP','bzHigh_iHP'});
            %  reshape([(Cells(num_cell).individual_trials_random(row_selected_ripples,1:length(edges_CCA))');  (Cells(num_cell).individual_trials(row_selected_ripples,1:length(edges_CCA))')],1,[]);
            Activitymatrix_temp(darab_neuron,:)= reshape([(Cells(num_cell).individual_trials_random(row_selected_ripples,1:length(edges_CCA))');  (Cells(num_cell).individual_trials(row_selected_ripples,1:length(edges_CCA))')],1,[]); %
            %reshape((Cells(num_cell).individual_trials(row_selected_ripples,1:length(edges_CCA))'),1,[]);%-mean(reshape(Cells(num_cell).individual_trials(row_selected_ripples,1:length(edges_CCA)-1)',1,[]));
            darab_neuron=darab_neuron+1;
        end


        Assembly_cellID_temp=[intersect(intersect(find(ismember({Cells.clusterNotes},{'single'})), find(ismember({Cells.Structre},structures_unique(structures_NUM)))),find(ismember([Cells.session],session_ID(nume_files))))];
        eval(['Assembly_cellID(nume_files).' char(structures_unique{structures_NUM}) '=Assembly_cellID_temp;'])

        eval(['Activitymatrix(nume_files).' char(structures_unique{structures_NUM}) '=Activitymatrix_temp;'])
        % Activitymatrix(nume_files).PrL=Activitymatrix_temp;

        eval(['AssemblyTemplates(nume_files).' char(structures_unique{structures_NUM}) '=assembly_patterns(Activitymatrix(nume_files).' char(structures_unique{structures_NUM}) ',opts);'])
        % AssemblyTemplates(nume_files).PrL=assembly_patterns(Activitymatrix(nume_files).PrL,opts);

        eval(['Activities(nume_files).' char(structures_unique{structures_NUM}) '=assembly_activity(AssemblyTemplates(nume_files).' char(structures_unique{structures_NUM}) ',Activitymatrix(nume_files).' char(structures_unique{structures_NUM}) ');' ])
        % Activities(nume_files).PrL = assembly_activity(AssemblyTemplates(nume_files).PrL,Activitymatrix(nume_files).PrL);

        for assembly_num=1:eval(['size(Activities(nume_files).' char(structures_unique{structures_NUM}) ',1)'])
            assembly_activity_temp(:,:,assembly_num)=eval(['reshape(Activities(nume_files).' char(structures_unique{structures_NUM}) '(assembly_num,:),len,sum(row_selected_ripples));']);
        end
        eval(['Assembly_activity(nume_files).' char(structures_unique{structures_NUM}) '=assembly_activity_temp']);

        clear Activitymatrix_temp assembly_activity Assembly_cellID_temp   assembly_activity_temp

    end
end

%% Plotting Assemblies

flag_plot=1;
color_lines=colormap('lines');

for nume_files=1:numel(path)
    disp(['Working on Session ID: ' num2str(path(nume_files).sessionID)])

    for structures_NUM=1:length(structures_unique)
        for assembly_num_total=1:eval(['(size(Assembly_activity(nume_files). ' char(structures_unique{structures_NUM}) ',3))'])
            trashold=eval(['2*mean(((std(abs(Assembly_activity(nume_files).' char(structures_unique{structures_NUM})  '(1:length(start_dtop_rip(1):start_dtop_rip(2)),:,assembly_num_total))))));']);
            eventSubset_temp(assembly_num_total,:)=eval(['(max(abs(Assembly_activity(nume_files).' char(structures_unique{structures_NUM}) '((start_dtop_rip(1):start_dtop_rip(2))+len/2,:,assembly_num_total)))>trashold);']);
            eventSubset_temp_before_binary(assembly_num_total,:)=eval(['(max(abs(Assembly_activity(nume_files).' char(structures_unique{structures_NUM}) '(1:length(start_dtop_rip(1):start_dtop_rip(2)),:,assembly_num_total)))>trashold);']);
            eventSubset_max(assembly_num_total,:)=eval(['(max(abs(Assembly_activity(nume_files).' char(structures_unique{structures_NUM}) '((start_dtop_rip(1):start_dtop_rip(2))+len/2,:,assembly_num_total))));']);
            eventSubset_max_before(assembly_num_total,:)=eval(['(max(abs(Assembly_activity(nume_files).' char(structures_unique{structures_NUM}) '(1:length(start_dtop_rip(1):start_dtop_rip(2)),:,assembly_num_total))));']);

        end

        temp_binary=sum(repmat((pow2(size(eventSubset_temp,1)-1:-1:0))',[1,size(eventSubset_temp,2)]).*eventSubset_temp);
        temp_binary_before=sum(repmat((pow2(size(eventSubset_temp_before_binary,1)-1:-1:0))',[1,size(eventSubset_temp_before_binary,2)]).*eventSubset_temp_before_binary);

        eval(['Assembly_decimal_ID(nume_files).' char(structures_unique{structures_NUM}) '= temp_binary;']);
        eval(['Assembly_decimal_rand_ID(nume_files).' char(structures_unique{structures_NUM}) '= temp_binary_before;']);

        eval(['Assembly_max_during_Ripple_ID(nume_files).' char(structures_unique{structures_NUM}) '= eventSubset_max;']);
        eval(['Assembly_max_random_ID(nume_files).' char(structures_unique{structures_NUM}) '=  eventSubset_max_before;']);

        if flag_plot==1
            fig_combination= figure;
            histogram(temp_binary_before(temp_binary_before>0),0.5:1:max(temp_binary)-0.5,'Normalization', 'probability');
            hold on;
            histogram(temp_binary(temp_binary>0),0.5:1:max(temp_binary)-0.5,'Normalization', 'probability');
            hold on;
        end

        trash=20;

        [N_peri,ind]= histc(temp_binary(temp_binary>0),1:1:max(temp_binary));
        edges_peri=(unique(ind));
        [N_pre,ind]= histc(temp_binary_before(temp_binary_before>0),1:1:max(temp_binary_before));
        edges_pre=(unique(ind));

        binary_num_peri=(find(N_peri>trash));
        binary_num_pre=(find(N_pre>trash));

        if flag_plot==1
            xticks(binary_num_peri)
            xticklabels({dec2bin(binary_num_peri)} )
            xtickangle(90)
            ylabel('Probability of significant Assembly Activition')
            xlabel('Pattern of Assembly activation before/during SWRs')
            title(char(structures_unique{structures_NUM}))

            fig=gcf;
            fig.PaperUnits = 'points';
            fig.Renderer='painters'
            fig.PaperPosition = [0 0 1200 800];
            fig.PaperSize = [1200 800];
            saveas(fig,[num2str(path(nume_files).sessionID) '_session_TBY35_Assembly_Probabilty_' char(structures_unique{structures_NUM})],'svg')
            saveas(fig,[num2str(path(nume_files).sessionID) '_session_TBY35_Assembly_Probabilty_' char(structures_unique{structures_NUM})],'tif')
            saveas(fig,[num2str(path(nume_files).sessionID) '_session_TBY35_Assembly_Probabilty_' char(structures_unique{structures_NUM})],'fig')
        end

        for i= 1:length(binary_num_peri)
            clear  Total_Ripple_selected_assembly   Ripples_selected_for_assembly_random Ripples_selected_for_assembly

            ripple_selected=find(ismember((Ripples(nume_files).ripples.eventID),{'bzHigh_dHP','bzLow_dHP','bzHigh_iHP','bzLow_iHP'}));
            eventSubset=intersect(ripple_selected,(find(temp_binary==binary_num_peri(i))));
            %             eventSubset_rand_temp=intersect(ripple_selected,(find(temp_binary~=binary_num_peri(i))));
            %             Ripples_selected_for_assembly=Ripples(nume_files).featureTable(eventSubset,1:60);
            %             entry = eventSubset_rand_temp(randperm(length(eventSubset_rand_temp)));
            %             eventSubset_rand=(entry(1:length(eventSubset)));
            %             Ripples_selected_for_assembly_random=Ripples(nume_files).featureTable(eventSubset_rand,1:60);
            %             Total_Ripple_selected_assembly=vertcat(Ripples_selected_for_assembly,Ripples_selected_for_assembly_random);
            %             Response=[repmat(1,[length(eventSubset) 1]); repmat(0,[length(eventSubset) 1])];




            if flag_plot==1
                figure(i)
                ax1=  subplot(1,3,1);
                for num_as=1:eval(['size(Assembly_activity(nume_files).' char(structures_unique{structures_NUM}) ',3)'])
                    %hold on; plot(edges_CCA,mean(assembly_activity_IL(:,eventSubset_IL,num_as),2),'-','LineWidth',2)
                    options.color_area = color_lines(num_as,:);% [243 169 114]./255;    % Orange theme
                    options.color_line = color_lines(num_as,:);%[236 112  22]./255;
                    temp_aerobar_random=eval(['Assembly_activity(nume_files).' char(structures_unique{structures_NUM}) '(1:length(edges_CCA),eventSubset,num_as);']);
                    plot_areaerrorbar(temp_aerobar_random',options);
                    temp_aerobar_random=[];
                    %      ylim([min(mean(Assembly_activity(nume_files).IL(1:41,eventSubset,num_as))), max(mean(Assembly_activity(nume_files).IL(1:41,eventSubset,num_as)))])
                    axis square
                end

                xlim([find(round(edges_CCA.*1000)==-500) find(round(edges_CCA.*1000)==500)])

                xticks([find(round(edges_CCA.*1000)==-500), find(round(edges_CCA.*1000)==-250), round(length(edges_CCA)/2),find(round(edges_CCA.*1000)==250), find(round(edges_CCA.*1000)==500)]);
                xticklabels([-500 -250,0,250,500])
                xlabel('Time [ms]')
                ylabel('Assembly expression strength')

                ax2= subplot(1,3,2);
                for num_as=1:eval(['size(Assembly_activity(nume_files).' char(structures_unique{structures_NUM}) ',3)'])
                    %hold on; plot(edges_CCA,mean(assembly_activity_IL(:,eventSubset_IL,num_as),2),'-','LineWidth',2)
                    options.color_area = color_lines(num_as,:);% [243 169 114]./255;    % Orange theme
                    options.color_line = color_lines(num_as,:);%[236 112  22]./255;
                    temp_aerobar_random=eval(['Assembly_activity(nume_files).' char(structures_unique{structures_NUM}) '(length(edges_CCA)+1:end,eventSubset,num_as);']);
                    plot_areaerrorbar(temp_aerobar_random',options);
                    temp_aerobar_random=[];


                    %ylim([-1 2])
                    axis square;
                end

                linkaxes([ax1,ax2],'y');

                xlim([find(round(edges_CCA.*1000)==-500) find(round(edges_CCA.*1000)==500)])
                xticks([find(round(edges_CCA.*1000)==-500), find(round(edges_CCA.*1000)==-250), round(length(edges_CCA)/2),find(round(edges_CCA.*1000)==250), find(round(edges_CCA.*1000)==500)]);
                xticklabels([-500 -250,0,250,500])
                xlabel('Time [ms]')
                %           ylabel('Assembly expression strength')
                %
                ax3=  subplot(1,3,3)

                %subplot(1,2,1)
                %mean(eventSubset_max)

                % figure;

                text_for_binary=dec2bin(binary_num_peri);

                violinplot((eventSubset_max(:,(temp_binary==binary_num_peri(i))))');
                title(['Assmebies pattern: ' (text_for_binary(i,:)) ' in ' char(structures_unique{structures_NUM})])
                ylabel('Assembly expression strength')
                xlabel('Assembly ID')
                % Value of the assembly strength ripple
                %  subplot(1,2,2)
                %  violinplot((eventSubset_max(:,find(temp_binary==binary_num(i-2))))')
                % title(['Number of ripples: ' num2str(numRipple(i-2))])
                % xlabel('Assembly ID')
                %  subplot(2,2,3)
                axis square
                fig=gcf;
                fig.PaperUnits = 'points';
                fig.Renderer='painters'
                fig.PaperPosition = [0 0 800 600];
                fig.PaperSize = [800 600];
                saveas(fig,[num2str(path(nume_files).sessionID) '_session_'  char(structures_unique{structures_NUM}) '_TBY37_Assembly_pattern_' num2str(dec2bin(binary_num_peri(i)))],'fig')
                saveas(fig,[num2str(path(nume_files).sessionID) '_session_'  char(structures_unique{structures_NUM}) '_TBY37_Assembly_pattern_' num2str(dec2bin(binary_num_peri(i)))],'svg')
                saveas(fig,[num2str(path(nume_files).sessionID) '_session_'  char(structures_unique{structures_NUM}) '_TBY37_Assembly_pattern_' text_for_binary(i,:)],'tif')

            end
            %saveas(fig,['TBY35_Assembly_patter_' num2str(dec2bin(binary_num_peri(i)))],'fig')

            close all
        end



        clear  eventSubset_temp  eventSubset_temp_before_binary   eventSubset_max  eventSubset_max_before
    end
end

%% OLD calculataion

% %% OLD calculataion
%
%
% handle=figure(2);
% handle.Renderer='opengl';
%
% subplot_counter=1
% for nume_files=1:numel(path)
%
%
%     for assembly_num_total=1:(size(Assembly_activity(nume_files).mPFC,3))
%         eventSubset_mPFC_temp(assembly_num_total,:)=(max(abs(Assembly_activity(nume_files).mPFC(start_dtop_rip(1):start_dtop_rip(2),:,assembly_num_total)))>5)';
%     end
%
%     ripple_selected=find(ismember((Ripples(nume_files).ripples.eventID),{'bzHigh_dHP','bzLow_dHP','bzLow_iHP','bzHigh_iHP'}))
%     eventSubset_mPFC=intersect(ripple_selected,(find(max(eventSubset_mPFC_temp))));
%
%     for assembly_num_total=1:(size(Assembly_activity(nume_files).RSC,3))
%         eventSubset_RSC_temp(assembly_num_total,:)=(max(abs(Assembly_activity(nume_files).RSC(start_dtop_rip(1):start_dtop_rip(2),:,assembly_num_total)))>5)';
%     end
%
%     if size(eventSubset_RSC_temp(assembly_num_total,:),1)>1
%         eventSubset_RSC=intersect(ripple_selected,(find(max(eventSubset_RSC_temp))));
%     else
%         eventSubset_RSC=intersect(ripple_selected,(find((eventSubset_RSC_temp))));
%     end
%
%     [eventSubset,ia,ib] =(setxor(eventSubset_mPFC, eventSubset_RSC));
%     eventSubset_mPFC=eventSubset_mPFC(ia);
%     eventSubset_RSC=eventSubset_RSC(ib);
%     intersect(eventSubset_mPFC,eventSubset_RSC)
%
%     EventSubset(nume_files).mPFC=eventSubset_mPFC;
%     EventSubset(nume_files).RSC=eventSubset_RSC;
%
%     subplot(5,2,subplot_counter)
%     for num_as=1:size(Assembly_activity(nume_files).mPFC,3)
%         %hold on; plot(edges_CCA,mean(assembly_activity_mPFC(:,eventSubset_mPFC,num_as),2),'-','LineWidth',2)
%         options.color_area = [243 169 114]./255;    % Orange theme
%         options.color_line = [236 112  22]./255;
%
%
%         plot_areaerrorbar((Assembly_activity(nume_files).mPFC(:,intersect(find(eventSubset_mPFC_temp(num_as,:)),eventSubset_mPFC'),num_as)'),options)
%         %ylim([min(mean(assembly_activity_mPFC(:,eventSubset_mPFC,num_as),2)), max(mean(assembly_activity_mPFC(:,eventSubset_mPFC,num_as),2))+0.1])
%     end
%
%     for num_as=1:size(Assembly_activity(nume_files).mPFC,3)
%         %   hold on; plot(edges_CCA,mean(assembly_activity_mPFC(:,eventSubset_RSC,num_as),2),'--','LineWidth',2)
%         options.color_area = [128 193 219]./255;    % Blue theme
%         options.color_line = [ 52 148 186]./255;
%         plot_areaerrorbar(Assembly_activity(nume_files).mPFC(:,eventSubset_RSC,num_as)',options)
%
%         %ylim([min(mean(assembly_activity_mPFC(:,eventSubset_mPFC,num_as),2)), max(mean(assembly_activity_mPFC(:,eventSubset_mPFC,num_as),2))+0.1])
%     end
%
%     xlim([find(round(edges_CCA.*1000)==-250) find(round(edges_CCA.*1000)==250)])
%     xticks([find(round(edges_CCA.*1000)==-250) round(length(edges_CCA)/2),find(round(edges_CCA.*1000)==250)]);
%     xticklabels([-250,0,250])
%
%     ylabel('assembly strength mPFC')
%     xlabel('time (msec)')
%     subtitle('Avaraged peri-SWR activation of mPFC assembly')
%     axis square
%
%     subplot_counter=subplot_counter+1
%
%     subplot(5,2,subplot_counter)
%     for num_as=1:size(Assembly_activity(nume_files).RSC,3)
%         %hold on; plot(edges_CCA,mean(assembly_activity_mPFC(:,eventSubset_mPFC,num_as),2),'-','LineWidth',2)
%         options.color_area = [243 169 114]./255;    % Orange theme
%         options.color_line = [236 112  22]./255;
%         plot_areaerrorbar(Assembly_activity(nume_files).RSC(:,eventSubset_mPFC,num_as)',options)
%         %ylim([min(mean(assembly_activity_mPFC(:,eventSubset_mPFC,num_as),2)), max(mean(assembly_activity_mPFC(:,eventSubset_mPFC,num_as),2))+0.1])
%     end
%
%     for num_as=1:size(Assembly_activity(nume_files).RSC,3)
%         %   hold on; plot(edges_CCA,mean(assembly_activity_mPFC(:,eventSubset_RSC,num_as),2),'--','LineWidth',2)
%         options.color_area = [128 193 219]./255;    % Blue theme
%         options.color_line = [ 52 148 186]./255;
%
%         plot_areaerrorbar((Assembly_activity(nume_files).RSC(:,intersect(find(eventSubset_RSC_temp(num_as,:)),eventSubset_RSC'),num_as)'),options)
%
%         %ylim([min(mean(assembly_activity_mPFC(:,eventSubset_mPFC,num_as),2)), max(mean(assembly_activity_mPFC(:,eventSubset_mPFC,num_as),2))+0.1])
%     end
%
%     xlim([find(round(edges_CCA.*1000)==-250) find(round(edges_CCA.*1000)==250)])
%     xticks([find(round(edges_CCA.*1000)==-250) round(length(edges_CCA)/2),find(round(edges_CCA.*1000)==250)]);
%     xticklabels([-250,0,250])
%
%     ylabel('assembly strength RSC')
%     xlabel('time msec)')
%     subtitle('Avaraged peri-SWR activation of RSC assembly')
%     axis square
%
%     subplot_counter=subplot_counter+1
%
%     clear eventSubset_mPFC_temp eventSubset_RSC_temp
% end
%
% fig=gcf;
% fig.PaperUnits = 'points';
% fig.Renderer='painters'
% fig.PaperPosition = [0 0 1500 1500];
% fig.PaperSize = [1500 1500];
% saveas(fig,['TBY37_Assembly'],'svg')
% saveas(fig,['TBY37_Assembly'],'fig')
load('Identified_Neurons.mat')
clear Identified_neruons
%%
% Identification
%this should give back the index from the Cells structure.

hol_session8=[Cells.session]==8;
hol_session9=[Cells.session]==9;
hol_session10=[Cells.session]==10;
hol_session11=[Cells.session]==11;
hol_session12=[Cells.session]==12;

ID_session=unique([Cells.session]);


for cells_num=130:size(Identified_Neurons,1)
    try
        if  Identified_Neurons.Session8(cells_num)==0
            Identified_neruons(cells_num,1)=0;
        elseif  Identified_Neurons.Session8(cells_num)~=0
            Identified_neruons(cells_num,1)=find([Cells(hol_session8).CellID_original]==Identified_Neurons.Session8(cells_num))+(find(hol_session8,1,'first')-1);
        end

        if  Identified_Neurons.Session9(cells_num)==0
            Identified_neruons(cells_num,2)=0;
        elseif  Identified_Neurons.Session9(cells_num)~=0
            Identified_neruons(cells_num,2)=find([Cells(hol_session9).CellID_original]==Identified_Neurons.Session9(cells_num))+(find(hol_session9,1,'first')-1);
        end

        if  Identified_Neurons.Session10(cells_num)==0
            Identified_neruons(cells_num,3)=0;
        elseif  Identified_Neurons.Session10(cells_num)~=0
            Identified_neruons(cells_num,3)=find([Cells(hol_session10).CellID_original]==Identified_Neurons.Session10(cells_num))+(find(hol_session10,1,'first')-1);
        end

        if  Identified_Neurons.Session11(cells_num)==0
            Identified_neruons(cells_num,4)=0;
        elseif  Identified_Neurons.Session11(cells_num)~=0
            Identified_neruons(cells_num,4)=find([Cells(hol_session11).CellID_original]==Identified_Neurons.Session11(cells_num))+(find(hol_session11,1,'first')-1);
        end

        if  Identified_Neurons.Session12(cells_num)==0
            Identified_neruons(cells_num,5)=0;
        elseif  Identified_Neurons.Session12(cells_num)~=0
            Identified_neruons(cells_num,5)=find([Cells(hol_session12).CellID_original]==Identified_Neurons.Session12(cells_num))+(find(hol_session12,1,'first')-1);
        end
    catch
        disp(['These cells are not exist in the Cells (Matlab):' num2str(cells_num) '+1 located in the unit_overview.ods'])
    end
end

%% STEM plotting for all structure

close all
for structures_NUM=1:length(structures_unique)
    figure(structures_NUM)

    %shift_subplot=[0 max((arrayfun(@(s)size(s.IL,1),Assembly_max_during_Ripple_ID))), max((arrayfun(@(s)size(s.IL,1),Assembly_max_during_Ripple_ID)))*2,max((arrayfun(@(s)size(s.IL,1),Assembly_max_during_Ripple_ID)))*3,max((arrayfun(@(s)size(s.IL,1),Assembly_max_during_Ripple_ID)))*4]
    shift_subplot=[0 eval(['max((arrayfun(@(s)size(s.' char(structures_unique(structures_NUM)) ',1),Assembly_max_during_Ripple_ID)))']), eval(['max((arrayfun(@(s)size(s.' char(structures_unique(structures_NUM)) ',1),Assembly_max_during_Ripple_ID)))*2']), eval(['max((arrayfun(@(s)size(s.' char(structures_unique(structures_NUM))  ',1),Assembly_max_during_Ripple_ID)))*3']),eval(['max((arrayfun(@(s)size(s.' char(structures_unique(structures_NUM)) ',1),Assembly_max_during_Ripple_ID)))*4'])];


    %max((arrayfun(@(s)size(s.IL,1),Assembly_max_during_Ripple_ID)))

    for nume_files=1:numel(path)
        s=1;

        %  cellId_text=[[(Cells(Assembly_cellID(nume_files).IL).CellID_original)];  [(Cells(Assembly_cellID(nume_files).IL).spikeSites)]; ];

        cellId_text=eval(['(Identified_neruons(ismember(Identified_neruons(:,nume_files),[(Assembly_cellID(nume_files).' char(structures_unique(structures_NUM)) ')]),:))']);
        cellId_text= cellId_text';

        [~,sorthely]=sort(cellId_text(nume_files,:))

        cellId_text=cellId_text(:,sorthely);

        temp_orig_ID=eval(['[(Assembly_cellID(nume_files).' char(structures_unique(structures_NUM)) ')]']);
       
        

try
        cellId_text(size(cellId_text,1)+1,:)= eval(['temp_orig_ID(ismember([(Assembly_cellID(nume_files).'  char(structures_unique(structures_NUM)) ')],Identified_neruons(:,nume_files)))']);
catch
    % cellId_text(size(cellId_text,1)+1,:)=temp_orig_ID;
end


        for assembly_num=1:eval(['size(AssemblyTemplates(nume_files).'  char(structures_unique(structures_NUM)) ',2)'])
            subplot_activity=subplot(5,eval(['max((arrayfun(@(s)size(s.'  char(structures_unique(structures_NUM)) ',1),Assembly_max_during_Ripple_ID)))']),s+shift_subplot(nume_files))

            significant=eval(['logical(sum([((AssemblyTemplates(nume_files).' char(structures_unique(structures_NUM)) '(:,assembly_num)))>1.5*std((AssemblyTemplates(nume_files).' char(structures_unique(structures_NUM)) '(:,assembly_num))), ((AssemblyTemplates(nume_files).' char(structures_unique(structures_NUM)) '(:,assembly_num)))<-1.5*std((AssemblyTemplates(nume_files).' char(structures_unique(structures_NUM)) '(:,assembly_num)))],2))']);
            % logical(sum([((AssemblyTemplates(nume_files).IL(:,assembly_num)))>1.5*std((AssemblyTemplates(nume_files).IL(:,assembly_num))), ((AssemblyTemplates(nume_files).IL(:,assembly_num)))<-1.5*std((AssemblyTemplates(nume_files).IL(:,assembly_num)))],2))


            [~,ia,~]=intersect(Identified_neruons(:,nume_files),eval(['Assembly_cellID(nume_files).' char(structures_unique(structures_NUM)) '(significant)']));
            eval(['Significant_Neurons_Assemblies(nume_files).' char(structures_unique(structures_NUM)) '(nume_files,assembly_num,:)={ia}']); %%%new
            eval(['Significant_Neurons_Assemblies(nume_files).' char(structures_unique(structures_NUM)) '_Original(nume_files,assembly_num,:)={Assembly_cellID(nume_files).' char(structures_unique(structures_NUM)) '(significant)}']);
            eval(['Significant_Neurons_Assemblies(nume_files).' char(structures_unique(structures_NUM)) '_significant(nume_files,assembly_num,:)=  {significant}']);
            clear ia;

            stem(find(significant),eval(['AssemblyTemplates(nume_files).' char(structures_unique(structures_NUM)) '(significant,assembly_num)']),'r-','filled')
            subtitle(num2str(assembly_num))
            hold on;
            stem(find(significant==0),eval(['AssemblyTemplates(nume_files).' char(structures_unique(structures_NUM)) '(~significant,assembly_num)']),'b-','filled')
            ylim([-1 1])
            xlim([0 length(cellId_text)+2])
            for num_id=1:size(cellId_text,2)
                %   cellId_Text(num_id)={(['N#' num2str(cellId_text(1,num_id)) ' S#' num2str(cellId_text(2,num_id)) ])}
                cellId_Text(num_id)={['S7 ID: ' num2str(cellId_text(1,num_id)) ' S8 ID: ' num2str(cellId_text(2,num_id)) ' S9 ID: ' num2str(cellId_text(3,num_id)) ' S10 ID: ' num2str(cellId_text(4,num_id)) ' S11 ID: ' num2str(cellId_text(5,num_id))   ' Cells ID ' num2str(cellId_text(end,num_id))]};
            end
            if assembly_num==1
                set(subplot_activity,'XTick',1:eval(['length([(Cells(Assembly_cellID(nume_files).' char(structures_unique(structures_NUM)) ').CellID_original)])']),'XTickLabel', cellId_Text,'FontSize',16');
            else
                %             set(subplot_activity,'XTick',[1:length([(Cells(Assembly_cellID(nume_files).IL).CellID_original)])],'XTickLabel',{(Cells(Assembly_cellID(nume_files).IL).CellID_original)},'FontSize',16')
                set(subplot_activity,'XTick',1:eval(['length([(Cells(Assembly_cellID(nume_files).' char(structures_unique(structures_NUM)) ').CellID_original)])']),'XTickLabel',   cellId_text(end,:),'FontSize',16')
            end
            view(90,-90)
            s=s+1;
        end
        clear   cellId_text  cellId_Text
    end

    sgtitle(['Assembly Template in '  char(structures_unique(structures_NUM))  ])

    fig=gcf
    fig.PaperUnits = 'points';
    fig.PaperPosition = [0 0 6000 4000];
    fig.PaperSize = [6000 4000];
    saveas(fig,['TBY37_Assembly_Template_' char(structures_unique(structures_NUM))],'svg')
    saveas(fig,['TBY37_Assembly_Template_' char(structures_unique(structures_NUM))],'fig')
end

close all

for structures_NUM=1:length(structures_unique)

    clear intersect_neurons_assembly_ID  intersect_neurons_assembly_percent  sorted_intresct_assemblies_similarity
    for darab_assembly=1:eval(['numel(Significant_Neurons_Assemblies(1).' char(structures_unique(structures_NUM)) '(1,:))'])
        % calculate the % of common cells within an assemblymatrix
        Ref_assembly=eval(['Significant_Neurons_Assemblies(1).' char(structures_unique(structures_NUM)) '{1,darab_assembly}(:,:)']);
        for session_n=1:numel(Significant_Neurons_Assemblies)
            for  darab_assembly_futo=1:size(eval(['Significant_Neurons_Assemblies(session_n).' char(structures_unique(structures_NUM))]),2)
                intersect_neurons_assembly_ID(session_n,darab_assembly_futo,darab_assembly)=eval(['{intersect(Ref_assembly,Significant_Neurons_Assemblies(session_n).' char(structures_unique(structures_NUM)) '{session_n, darab_assembly_futo}(:,:))}'])
                intersect_neurons_assembly_percent(session_n,darab_assembly_futo,darab_assembly)=eval(['length(intersect(Ref_assembly,Significant_Neurons_Assemblies(session_n).' char(structures_unique(structures_NUM)) '{session_n, darab_assembly_futo}(:,:)))./(length(Ref_assembly)+length(Significant_Neurons_Assemblies(session_n).' char(structures_unique(structures_NUM)) '{session_n, darab_assembly_futo}(:,:))-length(intersect(Ref_assembly,Significant_Neurons_Assemblies(session_n).' char(structures_unique(structures_NUM)) '{session_n, darab_assembly_futo}(:,:))))'])

            end
        end
    end
    intersect_neurons_assembly_percent(isnan(intersect_neurons_assembly_percent))=0;

    %Egy sorba helyezi a hasonlo halozatokat
%     for i=1:size(intersect_neurons_assembly_percent,3)
%         for session_n=1:size(intersect_neurons_assembly_percent,1)
%             sorted_intresct_assemblies_similarity(session_n,:,i)=sort(intersect_neurons_assembly_percent(session_n,:,i),'descend')
%         end
% 
%     end
     sorted_intresct_assemblies_similarity=intersect_neurons_assembly_percent;

    s_valtozo=([2:2:(size(sorted_intresct_assemblies_similarity,1)*2)])
    %cmap = colormap(lines);
    for assembly_num=1:size(sorted_intresct_assemblies_similarity,3)
        figure

        subplot((size(sorted_intresct_assemblies_similarity,1)),2,[1:2:(size(sorted_intresct_assemblies_similarity,1)*2)])
        heatmap(sorted_intresct_assemblies_similarity(:,:,assembly_num))

        title({'Similarity of the assemblies across sessions',['refAssembly ID:' num2str(assembly_num) ' Structure: '  char(structures_unique(structures_NUM))]})

        [~,max_sim_assembly_ID]=max(intersect_neurons_assembly_percent(:,:,assembly_num)')

        for nume_files=1:size(sorted_intresct_assemblies_similarity,1)
            subplot_activity=subplot((size(sorted_intresct_assemblies_similarity,1)),2, s_valtozo(nume_files))
            significant=eval(['[Significant_Neurons_Assemblies(nume_files).' char(structures_unique(structures_NUM)) '_significant{nume_files,max_sim_assembly_ID(nume_files)}(:,:)]'])
            if eval(['sum(AssemblyTemplates(nume_files).' char(structures_unique(structures_NUM))  '(significant,max_sim_assembly_ID(nume_files))<0)>1'])
                flip_szorzo=-1;
            else
                flip_szorzo=1;
            end
            stem(find(significant),eval(['[AssemblyTemplates(nume_files).' char(structures_unique(structures_NUM)) '(significant,max_sim_assembly_ID(nume_files))]']).*flip_szorzo,'r-','filled');
            set(subplot_activity,'XTick',find(significant),'XTickLabel',eval(['[Significant_Neurons_Assemblies(nume_files).' char(structures_unique(structures_NUM)) '{nume_files,max_sim_assembly_ID(nume_files)}(:,:)]']),'FontSize',16');
            hold on;
            stem(find(significant==0),eval(['[AssemblyTemplates(nume_files).' char(structures_unique(structures_NUM)) '(~significant,max_sim_assembly_ID(nume_files))]']).*flip_szorzo,'b-','filled')
            subtitle(num2str(max_sim_assembly_ID(nume_files)))
%             fig_ax=gca;
%               fig_ax.PlotBoxAspectRatio=([4,2,3]) ;
%           
%             view(90,-90)
        end
        fig=gcf
        fig.PaperUnits = 'points';
        fig.PaperPosition = [0 0 1400 600];
        fig.PaperSize = [1400 600];
        saveas(fig,['TBY37_Assembly_Tracking_'  char(structures_unique(structures_NUM)) '_assembly_' num2str(assembly_num)],'svg')
       saveas(fig,['TBY37_Assembly_Tracking_'  char(structures_unique(structures_NUM)) '_assembly_' num2str(assembly_num)],'tif')

    end
end


%
% %% percentage of matching for all, referenci is always the first session
% % Significant_Neurons_Assemblies.PrL;
%
%
% close all
% shift_subplot=[0 max((arrayfun(@(s)size(s.PrL,1),Assembly_max_during_Ripple_ID))), max((arrayfun(@(s)size(s.PrL,1),Assembly_max_during_Ripple_ID)))*2,max((arrayfun(@(s)size(s.PrL,1),Assembly_max_during_Ripple_ID)))*3,max((arrayfun(@(s)size(s.PrL,1),Assembly_max_during_Ripple_ID)))*4]
% %max((arrayfun(@(s)size(s.PrL,1),Assembly_max_during_Ripple_ID)))
%
% for nume_files=1:numel(path)
%     s=1;
%
%     %  cellId_text=[[(Cells(Assembly_cellID(nume_files).PrL).CellID_original)];  [(Cells(Assembly_cellID(nume_files).PrL).spikeSites)]; ];
%
%     cellId_text=(Identified_neruons(ismember(Identified_neruons(:,nume_files),[(Assembly_cellID(nume_files).PrL)]),:))';
%
%     [~,sorthely]=sort(cellId_text(nume_files,:))
%
%     cellId_text=cellId_text(:,sorthely);
%
%     temp_orig_ID=[(Assembly_cellID(nume_files).PrL)];
%     cellId_text(size(cellId_text,1)+1,:)= temp_orig_ID(ismember([(Assembly_cellID(nume_files).PrL)],Identified_neruons(:,nume_files)));
%
%     for assembly_num=1:size(AssemblyTemplates(nume_files).PrL,2)
%         subplot_activity=subplot(5,max((arrayfun(@(s)size(s.PrL,1),Assembly_max_during_Ripple_ID))),s+shift_subplot(nume_files))
%
%         significant=logical(sum([((AssemblyTemplates(nume_files).PrL(:,assembly_num)))>1.5*std((AssemblyTemplates(nume_files).PrL(:,assembly_num))), ((AssemblyTemplates(nume_files).PrL(:,assembly_num)))<-1.5*std((AssemblyTemplates(nume_files).PrL(:,assembly_num)))],2));
%
%         [~,ia,~]=intersect(Identified_neruons(:,nume_files),Assembly_cellID(nume_files).PrL(significant))
%         Significant_Neurons_Assemblies.PrL(nume_files,assembly_num,:)={ia}; %%%new
%         Significant_Neurons_Assemblies.PrL_Original(nume_files,assembly_num,:)={Assembly_cellID(nume_files).PrL(significant)};
%         Significant_Neurons_Assemblies.PrL_significant(nume_files,assembly_num,:)=  {significant};
%         clear ia;
%
%         stem(find(significant),AssemblyTemplates(nume_files).PrL(significant,assembly_num),'r-','filled')
%         subtitle(num2str(assembly_num))
%         hold on;
%         stem(find(significant==0),AssemblyTemplates(nume_files).PrL(~significant,assembly_num),'b-','filled')
%         ylim([-1 1])
%         xlim([0 length(cellId_text)+2])
%         for num_id=1:size(cellId_text,2)
%             %   cellId_Text(num_id)={(['N#' num2str(cellId_text(1,num_id)) ' S#' num2str(cellId_text(2,num_id)) ])}
%             cellId_Text(num_id)={['S7 ID: ' num2str(cellId_text(1,num_id)) ' S8 ID: ' num2str(cellId_text(2,num_id)) ' S9 ID: ' num2str(cellId_text(3,num_id)) ' S10 ID: ' num2str(cellId_text(4,num_id)) ' S11 ID: ' num2str(cellId_text(5,num_id))   ' Cells ID ' num2str(cellId_text(end,num_id))]};
%         end
%         if assembly_num==1
%             set(subplot_activity,'XTick',[1:length([(Cells(Assembly_cellID(nume_files).PrL).CellID_original)])],'XTickLabel', cellId_Text,'FontSize',16')
%         else
%             %             set(subplot_activity,'XTick',[1:length([(Cells(Assembly_cellID(nume_files).PrL).CellID_original)])],'XTickLabel',{(Cells(Assembly_cellID(nume_files).PrL).CellID_original)},'FontSize',16')
%             set(subplot_activity,'XTick',[1:length([(Cells(Assembly_cellID(nume_files).PrL).CellID_original)])],'XTickLabel',   cellId_text(end,:),'FontSize',16')
%         end
%         view(90,-90)
%         s=s+1;
%     end
%     clear   cellId_text  cellId_Text
% end
%
% sgtitle("Assembly Template in PrL")
%
% clear intersect_neurons_assembly_ID  intersect_neurons_assembly_percent
% for darab_assembly=1:numel(Significant_Neurons_Assemblies.PrL(1,:))
%     % calculate the % of common cells within an assemblymatrix
%     Ref_assembly=Significant_Neurons_Assemblies.PrL{1,darab_assembly}(:,:);
%     for session_n=1:size(Significant_Neurons_Assemblies.PrL(:,:),1)
%         for  darab_assembly_futo=1:size(Significant_Neurons_Assemblies.PrL,2)
%             intersect_neurons_assembly_ID(session_n,darab_assembly_futo,darab_assembly)={intersect(Ref_assembly,Significant_Neurons_Assemblies.PrL{session_n, darab_assembly_futo}(:,:))}
%             intersect_neurons_assembly_percent(session_n,darab_assembly_futo,darab_assembly)=length(intersect(Ref_assembly,Significant_Neurons_Assemblies.PrL{session_n, darab_assembly_futo}(:,:)))./(length(Ref_assembly)+length(Significant_Neurons_Assemblies.PrL{session_n, darab_assembly_futo}(:,:))-length(intersect(Ref_assembly,Significant_Neurons_Assemblies.PrL{session_n, darab_assembly_futo}(:,:))))
% %             intersect_assembly_ID(session_n,darab_assembly_futo)=darab_assembly_futo;
%         end
%     end
%
% end
%
% intersect_neurons_assembly_percent(isnan(intersect_neurons_assembly_percent))=0;
%
% %intersect_neurons_assembly_percent= intersect_neurons_assembly_percent.*100
%
% clear   sorted_intresct_assemblies_similarity
% for i=1:size(intersect_neurons_assembly_percent,3)
%     for session_n=1:size(intersect_neurons_assembly_percent,1)
%         sorted_intresct_assemblies_similarity(session_n,:,i)=sort(intersect_neurons_assembly_percent(session_n,:,i),'descend')
%     end
%
% end
%
%
%
% close all
% s_valtozo=([2:2:(size(sorted_intresct_assemblies_similarity,1)*2)])
% cmap = colormap(lines);
% for assembly_num=1:size(Assembly_max_during_Ripple_ID(1).PrL,1)
%     figure
%     subplot((size(sorted_intresct_assemblies_similarity,1)),2,[1:2:(size(sorted_intresct_assemblies_similarity,1)*2)])
%     heatmap(sorted_intresct_assemblies_similarity(:,:,assembly_num))
%     title(["Similarity of the assemblies across sessions","refAssembly ID:" num2str(assembly_num)])
%
%     [~,max_sim_assembly_ID]=max(intersect_neurons_assembly_percent(:,:,assembly_num)')
%
%     for nume_files=1:size(sorted_intresct_assemblies_similarity,1)
%         subplot_activity=subplot((size(sorted_intresct_assemblies_similarity,1)),2, s_valtozo(nume_files))
%
%
%         significant=[Significant_Neurons_Assemblies.PrL_significant{nume_files,max_sim_assembly_ID(nume_files)}(:,:)]
%
%
%         %         stem(find([Significant_Neurons_Assemblies.PrL_significant{nume_files,max_sim_assembly_ID(nume_files)}(:,:)]),AssemblyTemplates(nume_files).PrL([Significant_Neurons_Assemblies.PrL_significant{nume_files,max_sim_assembly_ID(nume_files)}(:,:)],s_valtozo(nume_files)),'Color',cmap( nume_files,:),'LineWidth',2,'MarkerFaceColor',cmap( nume_files,:))
%         %         hold on;
%         %         stem(find([Significant_Neurons_Assemblies.PrL_significant{nume_files,max_sim_assembly_ID(nume_files)}(:,:)]==0),AssemblyTemplates(nume_files).PrL(~[Significant_Neurons_Assemblies.PrL_significant{nume_files,max_sim_assembly_ID(nume_files)}(:,:)],s_valtozo(nume_files)),'b-','filled','LineWidth',2)
%         %   [~,~,holinter]=intersect(AssemblyTemplates(nume_files).PrL(significant,max_sim_assembly_ID(nume_files)),AssemblyTemplates(nume_files).PrL(:,max_sim_assembly_ID(nume_files)))
%
%
%         if sum(AssemblyTemplates(nume_files).PrL(significant,max_sim_assembly_ID(nume_files))<0)>1;
%             flip_szorzo=-1;
%         else
%             flip_szorzo=1;
%         end
%
%
%         stem(find(significant),[AssemblyTemplates(nume_files).PrL(significant,max_sim_assembly_ID(nume_files))].*flip_szorzo,'r-','filled')
%         set(subplot_activity,'XTick',find(significant),'XTickLabel',[Significant_Neurons_Assemblies.PrL{nume_files,max_sim_assembly_ID(nume_files)}(:,:)],'FontSize',16')
%
%         hold on;
%         stem(find(significant==0),[AssemblyTemplates(nume_files).PrL(~significant,max_sim_assembly_ID(nume_files))].*flip_szorzo,'b-','filled')
%         subtitle(num2str(max_sim_assembly_ID(nume_files)))
%
%
%
%     end
%
%     fig=gcf
%     fig.PaperUnits = 'points';
%     fig.PaperPosition = [0 0 1100 500];
%     fig.PaperSize = [1100 500];
%     saveas(fig,['TBY37_Assembly_Tracking_PrL_' num2str(assembly_num)],'fig')
%     saveas(fig,['TBY37_Assembly_Tracking_Prl_' num2str(assembly_num)],'tif')
%
% end
%
% %%
%
%
%
% figure(3)
%
% shift_subplot=[0 (size(AssemblyTemplates(1).mPFC,2)), (size(AssemblyTemplates(1).mPFC,2))*2,(size(AssemblyTemplates(1).mPFC,2))*3,(size(AssemblyTemplates(1).mPFC,2))*4]
% for nume_files=1:numel(path)
%     s=1;
%
%     %  cellId_text=[[(Cells(Assembly_cellID(nume_files).mPFC).CellID_original)];  [(Cells(Assembly_cellID(nume_files).mPFC).spikeSites)]; ];
%
%     cellId_text=(Identified_neruons(ismember(Identified_neruons(:,nume_files),[(Assembly_cellID(nume_files).mPFC)]),:))';
%
%     [~,sorthely]=sort(cellId_text(nume_files,:))
%
%     cellId_text=cellId_text(:,sorthely);
%
%     temp_orig_ID=[(Assembly_cellID(nume_files).mPFC)];
%     cellId_text(4,:)= temp_orig_ID(ismember([(Assembly_cellID(nume_files).mPFC)],Identified_neruons(:,nume_files)));
%
%     for assembly_num=1:size(AssemblyTemplates(nume_files).mPFC,2)
%         subplot_activity=subplot(5,(size(AssemblyTemplates(1).mPFC,2)),s+shift_subplot(nume_files))
%
%         significant=logical(sum([((AssemblyTemplates(nume_files).mPFC(:,assembly_num)))>1.5*std((AssemblyTemplates(nume_files).mPFC(:,assembly_num))), ((AssemblyTemplates(nume_files).mPFC(:,assembly_num)))<-1.5*std((AssemblyTemplates(nume_files).mPFC(:,assembly_num)))],2));
%         stem(find(significant),AssemblyTemplates(nume_files).mPFC(significant,assembly_num),'r-','filled')
%         hold on;
%         stem(find(significant==0),AssemblyTemplates(nume_files).mPFC(~significant,assembly_num),'b-','filled')
%         ylim([-1 1])
%         xlim([0 length(cellId_text)+2])
%         for num_id=1:size(cellId_text,2)
%
%
%             %   cellId_Text(num_id)={(['N#' num2str(cellId_text(1,num_id)) ' S#' num2str(cellId_text(2,num_id)) ])}
%             cellId_Text(num_id)={['S8 ID: ' num2str(cellId_text(1,num_id)) ' S10 ID: ' num2str(cellId_text(2,num_id)) ' S12 ID: ' num2str(cellId_text(3,num_id)),   ' Cells ID ' num2str(cellId_text(4,num_id))]};
%         end
%         if assembly_num==1
%             set(subplot_activity,'XTick',[1:length([(Cells(Assembly_cellID(nume_files).mPFC).CellID_original)])],'XTickLabel', cellId_Text,'FontSize',16')
%         else
%             %             set(subplot_activity,'XTick',[1:length([(Cells(Assembly_cellID(nume_files).mPFC).CellID_original)])],'XTickLabel',{(Cells(Assembly_cellID(nume_files).mPFC).CellID_original)},'FontSize',16')
%             set(subplot_activity,'XTick',[1:length([(Cells(Assembly_cellID(nume_files).mPFC).CellID_original)])],'XTickLabel',   cellId_text(4,:),'FontSize',16')
%         end
%         view(90,-90)
%         s=s+1;
%     end
%     clear   cellId_text  cellId_Text
% end
%
%
%
% fig=gcf
% fig.PaperUnits = 'points';
% fig.PaperPosition = [0 0 4000 2000];
% fig.PaperSize = [4000 2000];
% saveas(fig,['TBY37_Assembly_Template'],'svg')
% saveas(fig,['TBY37_Assembly_Template'],'fig')
%
%
% %%
%
% figure(4)
%
% shift_subplot=[0 (size(AssemblyTemplates(2).RSC,2)) ((size(AssemblyTemplates(2).RSC,2))*2)];
% for nume_files=1:numel(path)
%     s=1
%
%     cellId_text=(Identified_neruons(ismember(Identified_neruons(:,nume_files),[(Assembly_cellID(nume_files).RSC)]),:))';
%     [~,sorthely]=sort(cellId_text(nume_files,:))
%     cellId_text=cellId_text(:,sorthely);
%     temp_orig_ID=[(Assembly_cellID(nume_files).RSC)];
%     cellId_text(4,:)= temp_orig_ID(ismember([(Assembly_cellID(nume_files).RSC)],Identified_neruons(:,nume_files)));
%
%
%
%     for assembly_num=1:size(AssemblyTemplates(nume_files).RSC,2)
%         subplot_activity=subplot(3,(size(AssemblyTemplates(2).RSC,2)),s+shift_subplot(nume_files))
%
%         significant=logical(sum([((AssemblyTemplates(nume_files).RSC(:,assembly_num)))>1.5*std((AssemblyTemplates(nume_files).RSC(:,assembly_num))), ((AssemblyTemplates(nume_files).RSC(:,assembly_num)))<-1.5*std((AssemblyTemplates(nume_files).RSC(:,assembly_num)))],2));
%         stem(find(significant),AssemblyTemplates(nume_files).RSC(significant,assembly_num),'r-','filled')
%         hold on;
%         stem(find(significant==0),AssemblyTemplates(nume_files).RSC(~significant,assembly_num),'b-','filled')
%         ylim([-1 1])
%
%         for num_id=1:size(cellId_text,2)
%             cellId_Text(num_id)={['S8 ID: ' num2str(cellId_text(1,num_id)) ' S10 ID: ' num2str(cellId_text(2,num_id)) ' S12 ID: ' num2str(cellId_text(3,num_id)),   ' Cells ID ' num2str(cellId_text(4,num_id))]};
%         end
%         if assembly_num==1
%             set(subplot_activity,'XTick',[1:length([(Cells(Assembly_cellID(nume_files).RSC).CellID_original)])],'XTickLabel', cellId_Text,'FontSize',16')
%         else
%             set(subplot_activity,'XTick',[1:length([(Cells(Assembly_cellID(nume_files).RSC).CellID_original)])],'XTickLabel',   cellId_text(4,:),'FontSize',16')
%         end
%         view(90,-90)
%         s=s+1
%     end
%
% end
%
% fig=gcf
% fig.PaperUnits = 'points';
% fig.PaperPosition = [0 0 4000 2000];
% fig.PaperSize = [4000 2000];
% saveas(fig,['TBY37_Assembly_Template_RSC'],'svg')
% saveas(fig,['TBY37_Assembly_Template_RSC'],'fig')
%
% close all
%
%
%
% for nume_files=1:numel(path)
%     for assembly_num=1:size(AssemblyTemplates(nume_files).mPFC,2)
%         s=1;
%         figure(sscanf([num2str(nume_files) num2str(assembly_num)],'%i'))
%         significant=logical(sum([((AssemblyTemplates(nume_files).mPFC(:,assembly_num)))>1.5*std((AssemblyTemplates(nume_files).mPFC(:,assembly_num))), ((AssemblyTemplates(nume_files).mPFC(:,assembly_num)))<-1.5*std((AssemblyTemplates(nume_files).mPFC(:,assembly_num)))],2));
%         cellId_text=[[(Cells(Assembly_cellID(nume_files).mPFC).CellID_original)];  [(Cells(Assembly_cellID(nume_files).mPFC).spikeSites)]; Assembly_cellID(nume_files).mPFC];
%         for num_id=1:size(cellId_text,2)
%             cellId_Text(num_id)={(['N#' num2str(cellId_text(1,num_id)) ' S#' num2str(cellId_text(2,num_id)) ' C#' num2str(num2str(cellId_text(3,num_id)))])}
%         end
%
%         for darab_neuron=1:length(find(significant))
%             subplot_activity=subplot(1,length(find(significant)),s)
%             hely_szig=find(significant);
%             loc_wavform=Assembly_cellID(nume_files).mPFC
%             wavform=Cells(loc_wavform(1,hely_szig(darab_neuron))).wavform
%             time_wavform=((0:length(wavform)-1)/20)-1
%             plot(time_wavform,wavform,'LineWidth',3)
%             xlabel('Time [ms]')
%             ylabel('Amplitude [\muV]')
%             title(cellId_Text(1,hely_szig(darab_neuron)))
%             axis square
%             s=s+1
%         end
%
%
%         s=s+1
%
%
%         fig=gcf
%         fig.PaperUnits = 'points';
%         fig.PaperPosition = [0 0 800 800];
%         fig.PaperSize = [800 800];
%         saveas(fig,['TBY37_Wavforms_session_assembly_mPFC_' num2str(nume_files) num2str(assembly_num)],'svg')
%         saveas(fig,['TBY37_Wavforms_session_assembly_mPFC_' num2str(nume_files) num2str(assembly_num)],'fig')
%
%     end
%
% end
% close all
%
% for nume_files=1:numel(path)
%     for assembly_num=1:size(AssemblyTemplates(nume_files).RSC,2)
%         s=1;
%         figure(sscanf([num2str(nume_files) num2str(assembly_num)],'%i'))
%         significant=logical(sum([((AssemblyTemplates(nume_files).RSC(:,assembly_num)))>1.5*std((AssemblyTemplates(nume_files).RSC(:,assembly_num))), ((AssemblyTemplates(nume_files).RSC(:,assembly_num)))<-1.5*std((AssemblyTemplates(nume_files).RSC(:,assembly_num)))],2));
%         cellId_text=[[(Cells(Assembly_cellID(nume_files).RSC).CellID_original)];  [(Cells(Assembly_cellID(nume_files).RSC).spikeSites)]; Assembly_cellID(nume_files).RSC];
%         for num_id=1:size(cellId_text,2)
%             cellId_Text(num_id)={(['N#' num2str(cellId_text(1,num_id)) ' S#' num2str(cellId_text(2,num_id)) ' C#' num2str(num2str(cellId_text(3,num_id)))])}
%         end
%
%         for darab_neuron=1:length(find(significant))
%             subplot_activity=subplot(1,length(find(significant)),s)
%             hely_szig=find(significant);
%             loc_wavform=Assembly_cellID(nume_files).RSC
%             wavform=Cells(loc_wavform(1,hely_szig(darab_neuron))).wavform
%             time_wavform=((0:length(wavform)-1)/20)-1
%             plot(time_wavform,wavform,'LineWidth',3)
%             xlabel('Time [ms]')
%             ylabel('Amplitude [\muV]')
%             title(cellId_Text(1,hely_szig(darab_neuron)))
%             axis square
%             s=s+1
%         end
%
%
%         s=s+1
%
%
%         fig=gcf
%         fig.PaperUnits = 'points';
%         fig.PaperPosition = [0 0 800 800];
%         fig.PaperSize = [800 800];
%         saveas(fig,['TBY37_Wavforms_session_assembly_RSC_' num2str(nume_files) num2str(assembly_num)],'svg')
%         saveas(fig,['TBY37_Wavforms_session_assembly_RSC_' num2str(nume_files) num2str(assembly_num)],'fig')
%
%     end
%
% end
% close all
%
% figure(2)
% selected_ass_mPFC=[3,4,1];
% selected_ass_RSC=[1,1,1];
% subplot_counter=1;
% for nume_files=1:numel(path)
%
%
%     subplot(3,2,subplot_counter)
%     x1=mean(((Assembly_activity(nume_files).mPFC(start_dtop_rip(1)+1:start_dtop_rip(2)-1,EventSubset(nume_files).mPFC,selected_ass_mPFC(nume_files)))));
%     x2=mean(((Assembly_activity(nume_files).mPFC(start_dtop_rip(1)+1:start_dtop_rip(2)-1,EventSubset(nume_files).RSC,selected_ass_mPFC(nume_files)))));
%     [h,p,ci,stats]=ttest2(x1,x2)
%
%
%     xgroupdata=categorical([repmat(({'trig by mPFC spec.ripple'}),[length(x1) 1]); repmat(({'trig by RSC spec.ripple'}),[length(x2) 1])])
%     boxchart(flip(xgroupdata),[x2 x1]','GroupByColor',flip(xgroupdata),'Notch','on')
%     ylim([min([x1 x2])-0.2 max([x1 x2])-2])
%     ylabel('assembly strength/ripples')
%     title(['p-value: ' num2str(round(p,4))])
%
%     subplot_counter=subplot_counter+1;
%     subplot(3,2,subplot_counter)
%
%
%     x3=mean(((Assembly_activity(nume_files).RSC(start_dtop_rip(1)+1:start_dtop_rip(2)-1,EventSubset(nume_files).RSC,selected_ass_RSC(nume_files)))))
%     x4=mean(((Assembly_activity(nume_files).RSC(1:length(start_dtop_rip(1)+1:start_dtop_rip(2))-1,EventSubset(nume_files).mPFC,selected_ass_RSC(nume_files)))))
%     [h,p,ci,stats]=ttest2(x3,x4)
%
%     xgroupdata=categorical([repmat(({'trig by RSC spec.ripple'}),[length(x3) 1]); repmat(({'trig by mPFC spec.ripple'}),[length(x4) 1])])
%     boxchart(xgroupdata,[x3 x4],'GroupByColor',xgroupdata,'Notch','on')
%     ylim([min([x3 x4])-0.2 max([x3 x4])-2])
%     ylabel('assembly strength/ripples')
%     title(['p-value: ' num2str(round(p,4))])
%     subplot_counter=subplot_counter+1;
% end
%
% fig=figure(7)
% fig.PaperUnits = 'points';
% fig.PaperPosition = [0 0 1000 2000];
% fig.PaperSize = [1000 2000];
% fig.Renderer='painters'
% saveas(fig,['TBY34_summery'],'svg')
% saveas(fig,['TBY37_summery'],'fig')
%
%
%
%
% % figure
%
% % length(find(ismember(Ripples(1).ripples.eventID(EventSubset(1).mPFC),{'bzHigh_dHP','bzLow_dHP'})))
% % length(find(ismember(Ripples(1).ripples.eventID(EventSubset(1).mPFC),{'bzHigh_iHP','bzLow_iHP'})))
% %
% % length(find(ismember(Ripples(1).ripples.eventID(EventSubset(1).RSC),{'bzHigh_dHP','bzLow_dHP'})))
% % length(find(ismember(Ripples(1).ripples.eventID(EventSubset(1).RSC),{'bzHigh_iHP','bzLow_iHP'})))
%
%
%
