function [Data]=ImporterDAT(filename,TotalChannels,ChannelNumber)
% IMPORTERDAT Import ephys recordings from .dat files
%   [Data] = ImporterDat(filename,TotalChannels,ChannelNumber)
%   filename : location of input file
%   TotalChannels : number of channels in the input file
%   ChannelNumber : index of the channel to import

% maps existing file to memory and returns memory map. This allows access
% to files on disk as if they were in dynamic memory.
% data is accessible by calling outputObject.Data
a=memmapfile(filename, 'Format','int16');
% open file for binary read access and return integer file identifier
% here we create a new file by appending .new to the input file and open
% with write permission
% What is this used for ?
fid=fopen([filename '.new'], 'w');

Data=0;
block=1000000;

% create nonmodal dialoge box with waitbar and message
h = waitbar(0,'Import LFP from DAT file...');

% look up the size of the data stored in the recording file
total_per=size(a.data,1);

% iterate over start:step:stop_inclusive
% Start at ChannelNumber and add steps of size block * TotalChannels until
% you reach the size of the data set
for i=ChannelNumber:block*TotalChannels:size(a.data,1)

    % cap input chunk if we reach the end of file
    if (i+(block*TotalChannels) > size(a.data,1))
        b=a.data(i:TotalChannels:end);

    else
        % read in the data chunk from i to i+(blocksize *
        % TotalChannels)-1 in steps of TotalChannels
        % This reads the next block timesteps of ADC voltage samples
        % for channel ChannelNumber
        b=a.data(i:TotalChannels:i+(block* TotalChannels)-1);

    end
    % append to data
    Data=[Data; b];
    % update waitbar
    perc = round((size(Data,1)./(total_per./TotalChannels)),2)*100;
    waitbar(perc/100,h,sprintf('%.1f%% along...',perc))

end
Data(1)=[];
% clear variables and close waitbar
clear fid a perc
close(h)

end


