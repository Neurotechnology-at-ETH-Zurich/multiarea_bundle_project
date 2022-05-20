%this is a function to transafer clustered data for Cell Explorer
clc; 
clear;

%UI design for choosinfg bathpath
uiwait(msgbox('Please select your basepath.'))
basepath = uigetdir;
basepathPieces = regexp(basepath, filesep, 'split');
basename = basepathPieces{end};
ans_erase = questdlg('Erase previous data?');

if strcmp(ans_erase, 'Yes')
    %delete([basename,'.session.mat']);
    %delete([basename,'.spikes.cellinfo.mat']);
    delete(fullfile(basepath,[basename,'.cell_metrics.cellinfo.mat']));
    delete(fullfile(basepath,[basename,'.mono_res.cellinfo.mat']));
    delete(fullfile(basepath,[basename,'.session.mat']));
end

%% 1.1 loading data
%addpath(genpath("/Users/yuxuanhu/Projects/Spike_sort/rat/amplifier"));
%addpath(genpath("/Users/yuxuanhu/Documents/GitHub/CellExplorer"))
%basepath = '/Users/yuxuanhu/Projects/Spike_sort/rat/amplifier';

%% 1.3 creating a session 
session = sessionTemplate(basepath);

% And validate the required and optional fields
%validateSessionStruct(session);

%% 1.4 creating cell matrices
%Run the cell metrics pipeline 'ProcessCellMetrics' using the session struct as input
cell_metrics = ProcessCellMetrics('session', session);

% Visualize the cell metrics in CellExplorer
cell_metrics = CellExplorer('metrics',cell_metrics); 

