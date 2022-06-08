function varargout = show_spikes_info(varargin)
% SHOW_SPIKES_INFO MATLAB code for show_spikes_info.fig
%      SHOW_SPIKES_INFO, by itself, creates a new SHOW_SPIKES_INFO or raises the existing
%      singleton*.
%
%      H = SHOW_SPIKES_INFO returns the handle to a new SHOW_SPIKES_INFO or the handle to
%      the existing singleton*.
%
%      SHOW_SPIKES_INFO('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SHOW_SPIKES_INFO.M with the given input arguments.
%
%      SHOW_SPIKES_INFO('Property','Value',...) creates a new SHOW_SPIKES_INFO or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before show_spikes_info_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to show_spikes_info_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help show_spikes_info

% Last Modified by GUIDE v2.5 25-May-2022 17:17:17

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @show_spikes_info_OpeningFcn, ...
                   'gui_OutputFcn',  @show_spikes_info_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before show_spikes_info is made visible.
function show_spikes_info_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to show_spikes_info (see VARARGIN)

if ~exist('handles.cell_id','var')
uiwait(msgbox('Please select your basepath.'))
basepath = uigetdir;
basepathPieces = regexp(basepath, filesep, 'split');
basename = basepathPieces{end};

load([basepath,'/',basename,'.cell_metrics.cellinfo.mat']);
handles.cell_types = cell_metrics.putativeCellType;
handles.shank_id = cell_metrics.shankID;
handles.spikes = load([basepath,'/',basename,'.spikes.cellinfo.mat']);
handles.hmaps = load([basepath,'/',basename,'.hmaps.cellinfo.mat']);
handles.trace = load([basepath,'/',basename,'.trace.mat']);
handles.cell_id = 1;
end
handles.numcells = handles.spikes.spikes.numcells;

set(handles.id,'string',handles.cell_id);
handles.x = linspace(0,700,35);
handles.y = linspace(0,700,35);
show_heatmap(handles);
set(handles.cell_type,'string',handles.cell_types{handles.cell_id});
get_shank_region(handles);
% Choose default command line output for show_spikes_info
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes show_spikes_info wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = show_spikes_info_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in last.
function last_Callback(hObject, eventdata, handles)
% hObject    handle to last (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ids = handles.cell_id;
ids = ids-1;
if ids ==0
    ids = handles.numcells;
end
handles.cell_id = ids;
set(handles.id,'string',handles.cell_id);
show_heatmap(handles);
set(handles.cell_type,'string',handles.cell_types{handles.cell_id});
get_shank_region(handles);
guidata(hObject,handles);



% --- Executes on button press in next.
function next_Callback(hObject, eventdata, handles)
% hObject    handle to next (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ids = handles.cell_id;
ids = ids+1;
if (ids >78)
    ids = 1;
end
handles.cell_id = ids;
set(handles.id,'string',handles.cell_id);
show_heatmap(handles);
set(handles.cell_type,'string',handles.cell_types{handles.cell_id});
get_shank_region(handles);
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function id_CreateFcn(hObject, eventdata, handles)
% hObject    handle to id (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



function go_id_Callback(hObject, eventdata, handles)
% hObject    handle to go_id (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of go_id as text
%        str2double(get(hObject,'String')) returns contents of go_id as a double


% --- Executes during object creation, after setting all properties.
function go_id_CreateFcn(hObject, eventdata, handles)
% hObject    handle to go_id (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in go_to_id.
function go_to_id_Callback(hObject, eventdata, handles)
% hObject    handle to go_to_id (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ids = str2double(get(handles.go_id,'string'));
if (ids>0)&&(ids<78)
handles.cell_id = ids;
set(handles.id,'string',handles.cell_id);
show_heatmap(handles);
end
set(handles.cell_type,'string',handles.cell_types{handles.cell_id});
get_shank_region(handles);
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function cell_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cell_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function shankregion_CreateFcn(hObject, eventdata, handles)
% hObject    handle to shankregion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
