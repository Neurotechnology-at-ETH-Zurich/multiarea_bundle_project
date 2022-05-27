function [Data]=ImporterDAT_multi(filename,TotalChannels,channels)
% ImporterDat_multi Import multiple ephys recording channels from .dat files
%   filename        string  location of input file
%   TotalChannels   (1,)    number of channels in the input file
%   Channels        (:,1)   list of the channel numbers to import.

% maps existing file to memory and returns memory map. This allows access
% to files on disk as if they were in dynamic memory.
% data is accessible by calling outputObject.Data
a=memmapfile(filename, 'Format','int16');
% open file for binary read access and return integer file identifier
% here we create a new file by appending .new to the input file and open
% with write permission
% What is this used for ?
fid=fopen([filename '.new'], 'w');

chunk_size=1000000;

% create nonmodal dialoge box with waitbar and message
h = waitbar(0,'Import LFP from DAT file...');

% look up the size of the data stored in the recording file
total_per=size(a.data,1);

n_channels_import = numel(channels);
max_channel_import = channels(end);

% alocate read in array
b = zeros(n_channels_import,chunk_size);

% Initialize return variable
Data=zeros(n_channels_import,1);


% iterate over start:step:stop_inclusive
% Start at ChannelNumber and add steps of size block * TotalChannels until
% you reach the size of the data set
for chunk=0:chunk_size*TotalChannels:size(a.data,1)

    % cap input chunk if we reach the end of file
    if (chunk+max_channel_import+(chunk_size*TotalChannels) > size(a.data,1))
        % just grow b to the right size
        b = a.data(chunk+channels(1):TotalChannels:end);
        for index = 2:n_channels_import
            channel_index = channels(index);
            b  = [b  a.data(chunk+channel_index:TotalChannels:end)]; % add a row for each channel
            % TODO FIX LAST ARRAY CONCAT
        end
        b = b'; % transpose to channels x time

    else
        % read in the data chunk from i to i+(blocksize *
        % TotalChannels)-1 in steps of TotalChannels
        % This reads the next block timesteps of ADC voltage samples
        % for channel ChannelNumber
        for index = 1:n_channels_import
            channel_index = channels(index);
            b(index,:) = a.data(chunk+channel_index:TotalChannels:chunk+channel_index+(chunk_size* TotalChannels)-1);
        end

    end
    % append to data
    Data=[Data b];
    % update waitbar
    perc = round((size(Data,2)./(total_per./TotalChannels)),2)*100;
    waitbar(perc/100,h,sprintf('%.1f%% along...',perc))

end
% Remove first entry of zeros
Data(:,1)=[];
% clear variables and close waitbar
clear fid a perc
close(h)

end


