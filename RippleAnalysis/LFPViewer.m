function fig = LFPViewer(data,samplerate,ripple_timestamps,ripple_centers,ripple_classes,readonly)
%LFPViewer view local field potential recordings and annotated events
% data              (time,channels) recording data for visualization
% samplerate        (1,)    samplerate in Hz
% ripple_timestamps (:,2)   start and endpoints of ripples in seconds
% ripple_centers    (:,1)   center timepoints of ripples in seconds
% ripple_classes    (:,1)   categorical vector containing the class of each event.
% readonly          (1,)    boolean switch for enabling event editing operations


%Parent container
fig = uifigure;
fig.Name = 'LFP Viewer';

%App Data
appdata.displayInterval = [1 samplerate]; % start by displaying the first second
appdata.samplerate = samplerate; % Samplerate of the LFP signal
appdata.data = data;
appdata.framecount = size(data,1);
appdata.n_channels = size(data,2);
appdata.trace_spacing = 150; % vertical spacing of channels
appdata.ripple_timestamps = ripple_timestamps;
appdata.ripple_centers = ripple_centers;
assert(size(ripple_timestamps,1)==size(ripple_classes,1),'ripple_timestamps and ripple_classes do not agree in dimension');
assert(size(ripple_timestamps,1)==size(ripple_centers,1),'ripple_timestamps and ripple_centers do not agree in dimension');
assert(iscategorical(ripple_classes),'ripple_classes needs to be a categorical vector');
appdata.ripple_classes = ripple_classes; % category of ripple events
appdata.readonly = readonly;

%state flags and data fields for appending new ripples to the list
appdata.isSelectingStart = false;
appdata.selectedStartTimestamp = 0;
appdata.isSelectingEnd = false;
appdata.selectedEndTimestamp = 0;

%help message
if appdata.readonly
    appdata.helpMessage = 'READONLY | left/right arrow : prev/next frame, (shift) up/down arrow : prev/next (active) ripple, j : jump';
else
    appdata.helpMessage = strcat('left/right arrow : prev/next frame, (shift) up/down arrow : prev/next (active) ripple, j : jump',newline, ...
        ' a : add, d: delete, l : label, s : save, r/n : rename/new class');
end

%state flags and data fields for deleting existing ripples
appdata.isDeleting = false;
appdata.deletionTimestamp = 0;

% save appdata 
guidata(fig,appdata);

%layout
gl = uigridlayout(fig,[2 2]);
gl.RowHeight = {'fit','1x'};
gl.ColumnWidth = {'fit','1x'};

%Create components
fig_label = uilabel(gl);
fig_axes = uiaxes(gl);
fig_listbox = uilistbox(gl);

%Position Components
fig_label.Layout.Row = 1;
fig_label.Layout.Column = [1 2];
fig_axes.Layout.Row = 2;
fig_axes.Layout.Column = 2;
fig_listbox.Layout.Row = 2;
fig_listbox.Layout.Column = 1;

% Initialize GUI
fig_label.Text = appdata.helpMessage;
updateClasses(fig,fig_listbox); % initial call of updateClasses
drawLFP(fig,fig_axes); % trigger inital plot

%App Behaviour
fig.KeyPressFcn = {@processKeyPress, fig_axes, fig_label, fig_listbox};
fig.WindowButtonDownFcn = {@(src,event)uiresume(src)}; % call uiresume on click



end

function drawLFP(src, uiaxes)
% drawLFP trigger update of the lfp plot
% src       the parent uifigure
% uiaxes    the uiaxes to draw upon

% unpack guidata, extract data in display intervall and generate timestamps
appdata = guidata(src);
data_slice_start = appdata.displayInterval(1);
data_slice_stop = appdata.displayInterval(2);
data_slice = appdata.data(data_slice_start:data_slice_stop,:);
timestamps = (data_slice_start:data_slice_stop)/appdata.samplerate;
ripple_timestamps = appdata.ripple_timestamps;
ripple_centers = appdata.ripple_centers;

% verticaly space data by adding offset values
trace_offset = (1:appdata.n_channels) * appdata.trace_spacing;
onevector = ones(numel(timestamps),1);
offset =  (trace_offset(:) * onevector(:)' )'; % calculate outer product of vectors

% show the lfp traces in a stacked plot
plot(uiaxes,timestamps,double(data_slice)+offset,'Color',[0.5 0.5 0.5]);
xlabel(uiaxes,'time [s]');
xlim(uiaxes,[timestamps(1) timestamps(end)]);
ylabel(uiaxes,'');

% query y extent of plot
yc = uiaxes.YLim;

% visualize ripple intervals by rectangular annotations
for i = 1:size(ripple_timestamps,1)
    % get intervall and check visibility
    interval_timestamps = ripple_timestamps(i,:);
    interval_frames = interval_timestamps * appdata.samplerate;

    % test if either the interval start or end point is located inside the
    % time interval displayed in the plot
    if (interval_frames(1)>=data_slice_start && interval_frames(1)<= data_slice_stop) ...
            || (interval_frames(2)>=data_slice_start && interval_frames(2) <= data_slice_stop)

        % assemple patch rectangle coordinates
        xc = interval_timestamps;
        px = [ xc(1) xc(2) xc(2) xc(1)];
        py = [ yc(1) yc(1) yc(2) yc(2)];

        % extract class id and map to color
        class = appdata.ripple_classes(i);
        class_index = class == categories(appdata.ripple_classes); % get index of ripple class
        color = appdata.colormap(class_index,:);

        patch(uiaxes,px,py,color,'FaceAlpha',0.3,'EdgeColor','none');
        text(uiaxes,xc(1),yc(1)+200,class); % display class label in top left corner

    end
end

% calculate ripple center frames and check visibility
ripple_centers_frames = ripple_centers * appdata.samplerate;
visible_centers = (data_slice_start <= ripple_centers_frames) & (ripple_centers_frames <= data_slice_stop);
visible_centers_index = find(visible_centers);
% visualize ripple centers using vertical lines
for index = visible_centers_index' % iterate over entries of row vector
    % extract class id and map to color
    class = appdata.ripple_classes(index);
    class_index = (class == categories(appdata.ripple_classes)); % get index of ripple class
    color = appdata.colormap(class_index,:);
    xline(uiaxes,ripple_centers(index),'Color',color);
end

end

function processKeyPress(src,event,uiaxes,fig_label,fig_listbox)
% processKeyPress this callback executes whenever a key is pressed while
% LFPViewer is active. Here we process all user input, trigger the
% corresponding actions and maintain the program data fields.

% extract data from parent figure
appdata = guidata(src);
display_start = appdata.displayInterval(1);
display_end = appdata.displayInterval(2);
intervall_length = display_end - display_start;

switch event.Key
    case 'leftarrow' % shift frame back in time
        if display_start - intervall_length >= 0
            % shift by display intervall if possible
            display_start = display_start - intervall_length;
            display_end = display_end - intervall_length;
        else
            % otherwise just reset to start
            display_start = 1;
            display_end = intervall_length;
        end

    case 'rightarrow' % shift frame forward in time
        if display_end + intervall_length <= appdata.framecount
            display_start = display_start + intervall_length;
            display_end = display_end + intervall_length;
        else
            display_start = appdata.framecount - intervall_length;
            display_end = appdata.framecount;
        end
    case 'uparrow' % jump to the next (active) ripple event
        % calculate current center frame and timestamp
        centerFrame = (display_start+display_end)/2;
        centerTimestamp = centerFrame / appdata.samplerate;

        % find the next ripple by searching for the event with minimal
        % positive onset delay
        onset_delay = appdata.ripple_timestamps(:,1) - (centerTimestamp+0.001); % add small number to make delay negative for current ripple
        onset_delay(onset_delay <= 0) = inf; % assign infinite delay to all ripples preceding current timepoint
        
        % if modifier shift is active, restrict candidates to the currently
        % active class
        if ismember('shift',event.Modifier)
            active = appdata.ripple_classes == fig_listbox.Value;
            onset_delay(~active) = inf;
        end
        [next_ripple_delay, next_ripple_index] = min(onset_delay);

        if next_ripple_delay ~= inf
            % if it exists make it the new center
            future = appdata.ripple_timestamps(next_ripple_index,1) * appdata.samplerate; % get frame number of next ripple
            display_start = max([round(future - intervall_length/2), 1]); % restrict intervall to beginning if necessary
            display_end = min([display_start + intervall_length, appdata.framecount]); % restrict intervall to eof if we would extend past it
        end
    case 'downarrow' % jump to the previous (active) ripple event
        % calculate current center frame and timestamp
        centerFrame = (display_start+display_end)/2;
        centerTimestamp = centerFrame / appdata.samplerate;

        % find the previous ripple by searching for the event with least
        % negative onset delay.
        onset_delay = appdata.ripple_timestamps(:,1) - (centerTimestamp-0.001); % subtract small number to make delay positive for current ripple
        onset_delay(onset_delay >= 0) = -inf; % assign infinite delay to all ripples following current timepoint

        % if modifier shift is active, restrict candidates to the currently
        % active class
        if ismember('shift',event.Modifier)
            active = appdata.ripple_classes == fig_listbox.Value;
            onset_delay(~active) = -inf;
        end
        [next_ripple_delay, next_ripple_index] = max(onset_delay);

        if next_ripple_delay ~= -inf
            % if it exists make it the new center
            past = appdata.ripple_timestamps(next_ripple_index,1) * appdata.samplerate; % get frame number of next ripple
            display_start = max([round(past - intervall_length/2),1]); % restrict intervall to beginning if necessary
            display_end = min([display_start + intervall_length, appdata.framecount]); % restrict intervall to eof if we would extend past it
        end

    case 'a' % add ripple interval defined by start and end point
        % check if we have write permission
        if appdata.readonly
            return
        end
        % check that we are currently not selecting an end point
        if ~appdata.isSelectingEnd
            fig_label.Text = 'Click to select event start point';
            % set state flag for selecting start point and change mouse pointer
            appdata.isSelectingStart = true;
            src.Pointer = 'left';
            uiwait(src); % wait for click or key press
            appdata.selectedStartTimestamp = uiaxes.CurrentPoint(1,1);
            % reset help text and pointer appearance, set state flag
            src.Pointer = 'arrow';
            fig_label.Text = appdata.helpMessage;
            appdata.isSelectingStart = false;
            appdata.isSelectingEnd = true;
        else
            fig_label.Text = 'Click to select event end point';
            src.Pointer = 'left';
            uiwait(src); % wait for click and get timestamp
            appdata.selectedEndTimestamp = uiaxes.CurrentPoint(1,1);
            % reset app state and pointer appearance
            appdata.isSelectingEnd = false;
            src.Pointer = 'arrow';
            fig_label.Text = appdata.helpMessage;

            % append new ripple to all property lists of ripples
            new_interval = sort([appdata.selectedStartTimestamp appdata.selectedEndTimestamp]);
            appdata.ripple_timestamps(end+1,:) = new_interval;
            appdata.ripple_classes(end+1) = fig_listbox.Value; % newly created ripples are labeled with the currently active class
            appdata.ripple_centers(end+1) = sum(new_interval)/2;
        end
    case 'd' % delete all ripple intervals containing target point
        % check if we have write permission
        if appdata.readonly
            return
        end
        % get user input
        fig_label.Text = 'Click to delete interval containing point';
        src.Pointer = 'left';
        uiwait(src); % wait for click and get timestamp
        deletion_point = uiaxes.CurrentPoint(1,1);
        % reset app state and pointer appearance
        appdata.isSelectingEnd = false;
        src.Pointer = 'arrow';
        fig_label.Text = appdata.helpMessage;
        % find overlaping ripples and gather indices for deletion
        to_delete = [];
        for i = 1:size(appdata.ripple_timestamps,1)
            interval = appdata.ripple_timestamps(i,:);
            if deletion_point>=interval(1) && deletion_point<=interval(2)
                to_delete = [to_delete i]; % mark row for deletion
            end
        end
        % delete ripples from all property lists of ripples
        appdata.ripple_timestamps(to_delete,:) = []; % delete outside of for loop to prevent index errors
        appdata.ripple_classes(to_delete) = [];
        appdata.ripple_centers(to_delete) = [];
        % drop unused categories if any
        appdata.ripple_classes = removecats(appdata.ripple_classes);
    case 's' % save ripples to mat file
        % check if we have write permission
        if appdata.readonly
            return
        end
        ripple_timestamps = appdata.ripple_timestamps;
        ripple_classes = appdata.ripple_classes;
        [savefile,savepath] = uiputfile('curated_ripples.mat');
        if savefile ~= 0
            savepath = fullfile(savepath,savefile);
            save(savepath,"ripple_timestamps","ripple_classes");
            message = {'Saving Succesfull!',strcat('Saved ripple timestamps to ',savepath,'.mat')};
            uialert(src,message,'Saving Successfull','Icon','success');
        end
    case 'l' % label all rippples containing target point with active class
        % check if we have write permission
        if appdata.readonly
            return
        end
        classID = round(str2double(event.Key));
        % get user input
        fig_label.Text = strcat('Click to label intervall as',fig_listbox.Value);
        src.Pointer = 'left';
        uiwait(src); % wait for click and get timestamp
        selection_point = uiaxes.CurrentPoint(1,1);
        src.Pointer = 'arrow';
        fig_label.Text = appdata.helpMessage;
        % find overlaping ripples and gather indices for deletion
        to_label = [];
        for i = 1:size(appdata.ripple_timestamps,1)
            interval = appdata.ripple_timestamps(i,:);
            if selection_point>=interval(1) && selection_point<=interval(2)
                to_label = [to_label i]; % mark row for labeling
            end
        end
        % assign class ID to ripples
        appdata.ripple_classes(to_label) = fig_listbox.Value;

    case 'j' % jump to target time
        % open modal dialog to get user input
        prompt = 'Enter timepoint to jump to [s]';
        dialog_title = 'Jump Dialog';
        dims = [1 35];
        answer = inputdlg(prompt,dialog_title,dims);
        answer = str2double(answer);

        if answer >= 0 && answer <= (appdata.framecount / appdata.samplerate)
            % got valid timepoint
            targetTimepoint = answer;
            % make it the new center
            targetFrame = targetTimepoint * appdata.samplerate;
            display_start = max([round(targetFrame - intervall_length/2), 1]); % restrict intervall to beginning if necessary
            display_end = min([display_start + intervall_length, appdata.framecount]); % restrict intervall to eof if we would extend past it
        end
    case 'r' % rename active class
        % get active class
        active_class = fig_listbox.Value;
        % promt user to get new category label
        new_label = inputdlg('Enter new class label','Rename Class');
        if ~ismember(new_label,categories(appdata.ripple_classes))
            % rename category
            appdata.ripple_classes = renamecats(appdata.ripple_classes,active_class,new_label);
        else
            'LFPViewer rename : label allready in use'
        end
    case 'n' % create new class
        % promt user to get new category label
        new_label = inputdlg('Enter new class label','Add Class');
        if ~ismember(new_label,categories(appdata.ripple_classes))
            % rename category
            appdata.ripple_classes = addcats(appdata.ripple_classes,new_label);
        else
            'LFPViewer rename : label allready in use'
        end
end

% update application data and write back to parent
appdata.displayInterval = [display_start display_end];
guidata(src,appdata);

updateClasses(src,fig_listbox); % since this has side effects in appdata make sure that we write back modifications from above before invoking it

% update the plot
drawLFP(src,uiaxes);
end

function updateClasses(fig,fig_listbox)
% updateClasses update class category related information
% NOTE : this function has side effects in guidata !
appdata = guidata(fig);

% generate colormap for patches
appdata.classes_n = numel(categories(appdata.ripple_classes)); % count number of categories
appdata.colormap = hsv(appdata.classes_n); % (num_classes,3) list of rgb tripples

% update class list in listbox
fig_listbox.Items = categories(appdata.ripple_classes);
guidata(fig,appdata);
end
