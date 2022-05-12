classdef RippleCurator < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

    properties
        % user input data fields
        lfp_samplerate;
        lfp_signal;
        ripple_timestamps;
        ripple_centers;
        ripple_clusterID;
        ripple_feature_table;
        ripple_n;
        
        % GUI objects
        feature_window
        feature_instruction
        feature_Selector1
        feature_Selector2
        feature_clusterListbox
        feature_axes

        cluster_window
        cluster_axis

        % application data fields
        cluster_n
        cluster_labels
        cluster_colors
        rippleEpisodes

        default_message = 's : split | m : merge | d : delete | a : align | w : save';

    end

    methods
        function curator = RippleCurator(lfp_samplerate, lfp_signal, ripple_timestamps, ripple_centers, ripple_clusterID, ripple_feature_table)
            %UNTITLED Construct an instance of this class
            %   ripple_clusterID categorization of ripples into clusters.
            %   zero is reserved for false positives / deleted ripples.
            %   Cluster numbering must be contigous.

            % store input datafields
            curator.lfp_samplerate = lfp_samplerate;
            curator.lfp_signal = lfp_signal;
            curator.ripple_timestamps = ripple_timestamps;
            curator.ripple_centers = ripple_centers;
            curator.ripple_n = numel(ripple_centers);
            curator.ripple_clusterID = ripple_clusterID;
            curator.ripple_feature_table = ripple_feature_table;
           
            % Spin up the GUI
            curator.createGUI();

            % Initial evaluation of cluster related variables
            curator.updateClusters();

            % extract rippleEpisodes
            curator.extractRippleEpisodes();

            % Display welcome message
            curator.feature_instruction.Text = curator.default_message;
            
            % Populate Drop down menu with feature names
            curator.feature_Selector1.Items = curator.ripple_feature_table.Properties.VariableNames;
            curator.feature_Selector2.Items = curator.ripple_feature_table.Properties.VariableNames;
            
            % hook up list box state change callback
            curator.feature_clusterListbox.ValueChangedFcn = @curator.updateFeatureWindow;

            % hook up drop down state change callback
            curator.feature_Selector1.ValueChangedFcn = @curator.updateFeatureWindow;
            curator.feature_Selector2.ValueChangedFcn = @curator.updateFeatureWindow;

            % hook up key press callback for feature window
            curator.feature_window.KeyPressFcn = @curator.handleKeyPress;

            % Initial draw of feature window
            curator.updateFeatureWindow([],[]);
            % Initial draw of cluster window
            curator.updateClusterWindow([],[]);

        end

        function createGUI(curator)
            % CREATE FEATURE WINDOW
            curator.feature_window = uifigure;
            curator.feature_window.Name = 'Ripple Curator';
            
            % create layout
            layout = uigridlayout(curator.feature_window,[4,2]);
            layout.RowHeight = {'fit','fit','fit','1x'};
            layout.ColumnWidth = {'fit','1x'};

            % create ui components
            instruction = uilabel(layout);
            instruction.Layout.Row = 1;
            instruction.Layout.Column = [1 2];
            curator.feature_instruction = instruction;

            featureSelector1 = uidropdown(layout);
            featureSelector1.Layout.Row = 2;
            featureSelector1.Layout.Column = 2;
            curator.feature_Selector1 = featureSelector1;
            featureSelector2 = uidropdown(layout);
            featureSelector2.Layout.Row = 3;
            featureSelector2.Layout.Column = 2;
            curator.feature_Selector2 = featureSelector2;

            listbox = uilistbox(layout);
            listbox.Layout.Row = 4;
            listbox.Layout.Column = 1;
            listbox.Multiselect = 'on';
            curator.feature_clusterListbox = listbox;

            ax = uiaxes(layout);
            ax.Layout.Row = 4;
            ax.Layout.Column = 2;
            curator.feature_axes = ax;

            % CREATE CLUSTER WINDOW
            curator.cluster_window = uifigure;
            curator.cluster_window.Name = 'Ripple Clusters';

            layout = uigridlayout(curator.cluster_window,[1 1]);
            layout.RowHeight = {'1x'};
            layout.ColumnWidth = {'1x'};

            curator.cluster_axis = uiaxes(layout);
            curator.cluster_axis.Layout.Row = 1;
            curator.cluster_axis.Layout.Column = 1;

        end

        function updateFeatureWindow(curator, src, event)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            
            % retrieve clusters that are active / selected
            activeClusters = cellfun(@str2num,curator.feature_clusterListbox.Value);
%             display_colors = ones(curator.cluster_n,3)*0.7;
%             for cluster_index = activeClusters
%                 display_colors(cluster_index,:) = curator.cluster_colors(cluster_index,:);
%             end

            % retrieve features that should be visualized
            feature1name = curator.feature_Selector1.Value;
            feature2name = curator.feature_Selector2.Value;
            xcords = table2array(curator.ripple_feature_table(:,feature1name));
            ycords = table2array(curator.ripple_feature_table(:,feature2name));

            % exclude deleted ripples
            valid_ripples = curator.ripple_clusterID ~=0;
            % devide points into active and incative subsets
            active_ripples = ismember(curator.ripple_clusterID,activeClusters);
  
            % plot datapoints
            vi = valid_ripples & ~active_ripples;
            va = valid_ripples & active_ripples;
            plot(curator.feature_axes,0,0); % clear plot
            gscatter(curator.feature_axes,xcords(vi),ycords(vi),curator.ripple_clusterID(vi),[[0.7 0.7 0.7]],[],[])
            hold(curator.feature_axes,'on');
            gscatter(curator.feature_axes,xcords(va),ycords(va),curator.ripple_clusterID(va),curator.cluster_colors(activeClusters,:),[],[])
            hold(curator.feature_axes,'off');
            xlabel(curator.feature_axes,feature1name);
            ylabel(curator.feature_axes,feature2name);
            axis(curator.feature_axes,'padded')
        end

        function updateClusterWindow(curator, src, event)
            offset = 0;
            num_episodes_background = 25;
            eventWindowFrames = size(curator.rippleEpisodes,2);
            % clear figure
            plot(curator.cluster_axis,0,0)
            % iteratively add clusters to plot
            hold(curator.cluster_axis,'on');
            for cluster_index = 1:curator.cluster_n
                % collect cluster members
                subset = curator.ripple_clusterID == cluster_index;
                subsetEpisodes = curator.rippleEpisodes(subset,:);
                subsetMeanEpisode = squeeze(mean(subsetEpisodes,1));
                % increase offset by minimum value
                offset = offset + abs(min(subsetEpisodes,[],'all'));
                % plot background episodes (all or subsample)
                if sum(subset) <= num_episodes_background
                    plot(curator.cluster_axis,(1:eventWindowFrames)/curator.lfp_samplerate,...
                        subsetEpisodes+offset,'Color',[0.7 0.7 0.7]);
                else
                    subsubset = randsample(size(subsetEpisodes,1),num_episodes_background);
                    subsubsetEpisodes = subsetEpisodes(subsubset,:);
                    plot(curator.cluster_axis,(1:eventWindowFrames)/curator.lfp_samplerate,...
                        subsubsetEpisodes+offset,'Color',[0.7 0.7 0.7]);
                end
                % plot mean episode
                plot(curator.cluster_axis,(1:eventWindowFrames)/curator.lfp_samplerate,...
                    subsetMeanEpisode+offset,'Color',curator.cluster_colors(cluster_index,:));
                
                % increase offset by maximum value
                offset = offset + abs(max(subsetEpisodes,[],'all'));
            end
            hold(curator.cluster_axis,'off');
            axis(curator.cluster_axis,'padded');
            xlabel(curator.cluster_axis,'time (s)');
            ylabel(curator.cluster_axis,'LFP Signal');
        end

        function handleKeyPress(curator, src, event)
            % handle key press events in the feature window
            switch event.Key
                case 's' % SELECTION
                    curator.feature_instruction.Text = 'Select points to create new cluster';
                    selection = curator.makeSelection(); % select a subset of the active points
                    % create a new cluster from the selection
                    curator.ripple_clusterID(selection) = curator.cluster_n+1;
                    curator.updateClusters();

                case 'd' % DELETE
                    curator.deleteClusters();
                    curator.updateClusters();
                case 'a' % ALIGN
                    curator.alignClusters();
                case 'm' % MERGE
                    curator.mergeClusters();
                    curator.updateClusters();
                case 'w' % SAVE
                    ripple_timestamps = curator.ripple_timestamps;
                    ripple_center_timestampls = curator.ripple_centers;
                    ripple_classes = curator.ripple_clusterID;
                    savepath = uiputfile('curated_ripples.mat');
                    if savepath ~= 0
                        save(savepath,"ripple_timestamps","ripple_classes","ripple_center_timestampls");
                        message = {'Saving Succesfull!',strcat('Saved ripple timestamps to ',savepath,'.mat')};
                        uialert(src,message,'Saving Successfull','Icon','success');
                    end
                    
            end
            % trigger update of application windows
            curator.feature_instruction.Text = curator.default_message;
            curator.updateFeatureWindow();
            curator.updateClusterWindow();
        end

        function extractRippleEpisodes(curator)
            % extract fixed length ripple episodes centered at user
            % specified timepoints
            % Set Window Length for Event Alignment
            eventWindow = 0.15; % [s] duration of event window
            eventWindowFrames = eventWindow*curator.lfp_samplerate; % convert to frame number
            eventWindowHalfFrames = round(eventWindowFrames /2); % extend central frame by this amount into both directions
            eventWindowFrames = 2*eventWindowHalfFrames + 1; % effective episode length
            % extract ripple episodes 
            curator.rippleEpisodes = zeros(curator.ripple_n,eventWindowFrames);
            % extract signals 
            for i = 1:curator.ripple_n
                center_frame = round(curator.ripple_centers(i)*curator.lfp_samplerate);
                startFrame = center_frame - eventWindowHalfFrames;
                endFrame = center_frame + eventWindowHalfFrames;
                curator.rippleEpisodes(i,:) = curator.lfp_signal(startFrame:endFrame);
            end
        end

        function updateClusters(curator)
            % update cluster related variables
            % calculate class number and assign colors
            curator.cluster_n = numel(setdiff(unique(curator.ripple_clusterID), [ 0 ] )); % do not count label zero as a class
            curator.cluster_colors = hsv(curator.cluster_n);
            % create class labels from integers
            curator.cluster_labels = {};
            for cluster_index = 1:curator.cluster_n
                curator.cluster_labels{cluster_index} = num2str(cluster_index);
            end

            % Populate List box with class labels
            curator.feature_clusterListbox.Items = curator.cluster_labels;
            curator.feature_clusterListbox.Value = curator.cluster_labels; % select all by default
        end

        function alignClusters(curator)
            % Refine centerpoints within each cluster
            % Set Window Length for Event Alignment
            eventWindow = 0.2; % [s] duration of event window
            eventWindowFrames = eventWindow*curator.lfp_samplerate; % convert to frame number
            eventWindowHalfFrames = round(eventWindowFrames /2); % extend central frame by this amount into both directions
            eventWindowFrames = 2*eventWindowHalfFrames + 1; % effective episode length
            % iterate over clusters
            for cluster_index = 1:curator.cluster_n
                members = curator.ripple_clusterID == cluster_index;
                center_frames = round(curator.ripple_centers(members)*curator.lfp_samplerate);
                % use cross correlation based alignment to line up the
                % ripple episodes of each cluster
                updatedCenters = alignRipples(curator.lfp_signal,center_frames,eventWindowFrames);
                curator.ripple_centers(members) = updatedCenters/curator.lfp_samplerate;
            end
            % extract ripple episodes centered at new timepoints
            curator.extractRippleEpisodes();
        end

        function selectedPoints = makeSelection(curator)
            % select a subset of the active points in feature window
            active_clusters = cellfun(@str2num,curator.feature_clusterListbox.Value);
            active = ismember(curator.ripple_clusterID,active_clusters);
            % retrieve coordinate that are displayed in feature view
            feature1name = curator.feature_Selector1.Value;
            feature2name = curator.feature_Selector2.Value;
            xcords = table2array(curator.ripple_feature_table(:,feature1name));
            ycords = table2array(curator.ripple_feature_table(:,feature2name));
            % make a selection
            selector = drawpolygon(curator.feature_axes);
            selection = inpolygon(xcords,ycords,selector.Position(:,1),selector.Position(:,2));
            % return points that are both active and selected
            selectedPoints = active & selection;
        end

        function deleteClusters(curator)
            % get active clusters
            active_clusters = cellfun(@str2num,curator.feature_clusterListbox.Value);
            % promt user for confirmation
            message = strcat('Delete clusters ',num2str(active_clusters),' ?');
            selection = uiconfirm(curator.feature_window,message,'Confirm Delete');
            if strcmp(selection,'OK')
                % delete selected clusters
                members = ismember(curator.ripple_clusterID,active_clusters);
                curator.ripple_clusterID(members) = 0; % assign to delete cluster
                % close gaps in numbering scheme by decreasing remaining
                % cluster indices greater than the deleted ones
                shift = zeros(curator.ripple_n,1);
                for deletedCluster = active_clusters
                    shift = shift + (curator.ripple_clusterID>deletedCluster);
                end
                curator.ripple_clusterID = curator.ripple_clusterID - shift;
            end
        end

        function mergeClusters(curator)
            % get active clusters
            active_clusters = cellfun(@str2num,curator.feature_clusterListbox.Value);
            % check that more than one cluster is selected
            if numel(active_clusters)>1
                % promt user for confirmation
                message = strcat('Merge clusters ',num2str(active_clusters),' ?');
                selection = uiconfirm(curator.feature_window,message,'Confirm Delete');
                if strcmp(selection,'OK')
                    % merge selected cluster
                    target_cluster = min(active_clusters); % merge into the smallest cluster index
                    move_clusters = setdiff(active_clusters,target_cluster);

                    members = ismember(curator.ripple_clusterID,move_clusters);
                    curator.ripple_clusterID(members) = target_cluster; % assign to target cluster
                    % close gaps in numbering scheme by decreasing remaining
                    % cluster indices greater than the moved ones
                    shift = zeros(curator.ripple_n,1);
                    for movedCluster = move_clusters
                        shift = shift + (curator.ripple_clusterID>movedCluster);
                    end
                    curator.ripple_clusterID = curator.ripple_clusterID - shift;
                end
            else
                'curator.mergeClusters : select more then one cluster'
            end
        end
    end
end