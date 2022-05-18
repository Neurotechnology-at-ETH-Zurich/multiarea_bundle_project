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

        lfp_viewer

        % application data fields
        cluster_n % the number of clusters including deleted ripples
        cluster_labels % cell array of strings holding cluster category names
        cluster_colors % mapping of cluster_label -> rgb tripplet
        rippleEpisodes

        default_message = 's : split | m : merge | d : delete | a : align | r : rename | w : save';

    end

    methods
        function curator = RippleCurator(lfp_samplerate, lfp_signal, ripple_timestamps, ripple_centers, ripple_categories, ripple_feature_table)
            %UNTITLED Construct an instance of this class
            %   ripple_clusterID (:,1) categorical array of cluster
            %   assignments. All ripples deleted in the session will be
            %   assigned to category 'deleted'.
            %   zero is reserved for false positives / deleted ripples.
            %   Cluster numbering must be contigous.
            
            % basic compatibilty checks
            assert(iscategorical(ripple_categories))

            % store input datafields
            curator.lfp_samplerate = lfp_samplerate;
            curator.lfp_signal = lfp_signal;
            curator.ripple_timestamps = ripple_timestamps;
            curator.ripple_centers = ripple_centers;
            curator.ripple_n = numel(ripple_centers);
            curator.ripple_clusterID = ripple_categories;
            curator.ripple_feature_table = ripple_feature_table;
           
            % Spin up the GUI
            curator.createGUI();

            % Initial evaluation of cluster related variables
            curator.ripple_clusterID = addcats(curator.ripple_clusterID,'deleted');
            curator.updateClusters();

            % start up the lfp viewer in read only mode
            curator.lfp_viewer = LFPViewer(lfp_signal,lfp_samplerate,ripple_timestamps,ripple_centers,ripple_categories,"",true);

            % Arrange windows on screen
            screen_size = get(0,'ScreenSize');
            width = screen_size(3);
            height = screen_size(4);
            curator.feature_window.Position = [1 height*0.3 width*0.6 height*0.65];
            curator.cluster_window.Position = [width*0.65 height*0.05 width*0.3 height*0.9];
            curator.lfp_viewer.Position = [1 height*0.05 width*0.6 height*0.2];
            
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

            % hook up callback for close request
            curator.feature_window.CloseRequestFcn = @curator.exitSession;
            curator.cluster_window.CloseRequestFcn = @curator.exitSession;

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
            activeClusters = curator.feature_clusterListbox.Value;
            
            % retrieve features that should be visualized
            feature1name = curator.feature_Selector1.Value;
            feature2name = curator.feature_Selector2.Value;
            xcords = table2array(curator.ripple_feature_table(:,feature1name));
            ycords = table2array(curator.ripple_feature_table(:,feature2name));

            % exclude deleted ripples
            valid_ripples = curator.ripple_clusterID ~= 'deleted';
            % devide points into active and incative subsets
            active_ripples = ismember(curator.ripple_clusterID,activeClusters);
  
            % plot datapoints
            vi = valid_ripples & ~active_ripples;
            va = valid_ripples & active_ripples;
            plot(curator.feature_axes,0,0); % clear plot
            gscatter(curator.feature_axes,xcords(vi),ycords(vi),curator.ripple_clusterID(vi),[[0.7 0.7 0.7]],[],[])
            hold(curator.feature_axes,'on');
            gscatter(curator.feature_axes,xcords(va),ycords(va),curator.ripple_clusterID(va));
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
            
            for cluster = setdiff(categories(curator.ripple_clusterID),{'deleted'},'stable')'
                % collect cluster members
                subset = curator.ripple_clusterID == cluster;
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
                % get cluster index
                cluster_index = strcmp(cluster,curator.cluster_labels);
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
                    curator.splitClusters();
                    curator.updateClusters();
                    curator.updateLFPViewer();
                case 'd' % DELETE
                    curator.deleteClusters();
                    curator.updateClusters();
                    curator.updateLFPViewer();
                case 'a' % ALIGN
                    curator.alignClusters();
                    curator.updateLFPViewer();
                case 'm' % MERGE
                    curator.mergeClusters();
                    curator.updateClusters();
                    curator.updateLFPViewer();
                case 'r' % RENAME
                    curator.renameCluster();
                    curator.updateClusters();
                    curator.updateLFPViewer();
                case 'w' % SAVE
                    ripple_timestamps = curator.ripple_timestamps;
                    ripple_center_timestamps = curator.ripple_centers;
                    ripple_classes = curator.ripple_clusterID;
                    savepath = uiputfile('curated_ripples.mat');
                    if savepath ~= 0
                        save(savepath,"ripple_timestamps","ripple_classes","ripple_center_timestamps");
                        message = {'Saving Succesfull!',strcat('Saved ripple timestamps to ',savepath,'.mat')};
                        uialert(src,message,'Saving Successfull','Icon','success');
                    end
                    
            end
            % trigger update of application windows
            curator.feature_instruction.Text = curator.default_message;
            curator.updateFeatureWindow();
            curator.updateClusterWindow();
        end

        function exitSession(curator, src, event)
            % close application windows
            selection = uiconfirm(src,'Exit Session ?','Confirm Close');
            switch selection
                case 'OK'
                    delete(curator.cluster_window);
                    delete(curator.feature_window);
                    delete(curator.lfp_viewer);
                case 'Cancel'
            end
            % Send the most relevant data to the workspace as a session
            % struct
            session.ripple_timestamps = curator.ripple_timestamps;
            session.ripple_center_timestamps = curator.ripple_centers;
            session.ripple_clusters = curator.ripple_clusterID;
            assignin("base","RippleCuratorSession",session);
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
        
        function splitClusters(curator)
            % TODO implement spli / merge / delete operations using
            % categorical class variables
            curator.feature_instruction.Text = 'Select points to create new cluster';
            selection = curator.makeSelection(); % select a subset of the active points
            % Generate a new class label that is not yet occupied
            i = 1;
            new_cluster_label = strcat('selection_',num2str(i));
            while ismember(new_cluster_label,curator.cluster_labels)
                i = i + 1;
                new_cluster_label = strcat('selection_',num2str(i));
            end
            % create a new cluster from the selection
            curator.ripple_clusterID(selection) = new_cluster_label;
        end

        function updateClusters(curator)
            % update cluster related variables
            curator.cluster_labels = categories(curator.ripple_clusterID);
            curator.cluster_n = numel(curator.cluster_labels);
            curator.cluster_colors = hsv(curator.cluster_n);
           
            % Populate List box with class labels
            curator.feature_clusterListbox.Items = setdiff(curator.cluster_labels,'deleted','stable');
            curator.feature_clusterListbox.Value = setdiff(curator.cluster_labels,'deleted','stable'); % select all by default
        end

        function updateLFPViewer(curator)
            % extract appdata
            data = guidata(curator.lfp_viewer);
            % update cluster related data fields
            data.ripple_timestamps = curator.ripple_timestamps;
            data.ripple_centers = curator.ripple_centers;
            data.ripple_classes = curator.ripple_clusterID;
            % write back to appdata
            guidata(curator.lfp_viewer,data);
        end

        function alignClusters(curator)
            % Refine centerpoints within each cluster
            % Set Window Length for Event Alignment
            eventWindow = 0.2; % [s] duration of event window
            eventWindowFrames = eventWindow*curator.lfp_samplerate; % convert to frame number
            eventWindowHalfFrames = round(eventWindowFrames /2); % extend central frame by this amount into both directions
            eventWindowFrames = 2*eventWindowHalfFrames + 1; % effective episode length
            % iterate over clusters
            for cluster = setdiff(categories(curator.ripple_clusterID),{'deleted'},'stable')'
                members = curator.ripple_clusterID == cluster;
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
            active_clusters = curator.feature_clusterListbox.Value;
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
            active_clusters = curator.feature_clusterListbox.Value;
            % promt user for confirmation
            message = strcat('Delete clusters?  ',strjoin(active_clusters));
            selection = uiconfirm(curator.feature_window,message,'Confirm Delete');
            if strcmp(selection,'OK')
                % delete selected clusters
                members = ismember(curator.ripple_clusterID,active_clusters);
                curator.ripple_clusterID(members) = 'deleted'; % assign to delete cluster
                % remove unused categories
                curator.ripple_clusterID = removecats(curator.ripple_clusterID);
            end
        end

        function mergeClusters(curator)
            % get active clusters
            active_clusters = curator.feature_clusterListbox.Value;
            % check that more than one cluster is selected
            if numel(active_clusters)>1
                % promt user for confirmation
                message = strcat('Merge clusters?  ',strjoin(active_clusters));
                selection = uiconfirm(curator.feature_window,message,'Confirm Merge');
                if strcmp(selection,'OK')
                    % merge selected cluster
                    target_cluster = active_clusters{1}; % merge into the first active category
                    move_clusters = setdiff(active_clusters,target_cluster);

                    members = ismember(curator.ripple_clusterID,move_clusters);
                    curator.ripple_clusterID(members) = target_cluster; % assign to target cluster
                    % remove unused categories
                    curator.ripple_clusterID = removecats(curator.ripple_clusterID);
                end
            else
                'curator.mergeClusters : select more then one cluster'
            end
        end

        function renameCluster(curator)
            % get active clusters
            active_clusters = curator.feature_clusterListbox.Value;
            % check that only one cluster is selected
            if numel(active_clusters) == 1
                % promt user to get new category label
                new_label = inputdlg('Enter new cluster label','Rename Cluster');
                if ~ismember(new_label,curator.cluster_labels)
                    % rename category
                    curator.ripple_clusterID = renamecats(curator.ripple_clusterID,active_clusters,new_label);
                else
                    'curator.renameCluster : label allready in use'
                end
            else
                'curator.renameCluster : select only one cluster'
            end
        end
    end
end