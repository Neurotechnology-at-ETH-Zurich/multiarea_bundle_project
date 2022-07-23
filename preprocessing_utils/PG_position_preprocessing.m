%% Import Coordinate from DLC csv
%written by Peter Gombkoto Ph.D ETH - Neurotechnology Group - pgombkoto@ethz.ch

% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 7);

% Specify range and delimiter
opts.DataLines = [4, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["frame", "RostralX", "RostralY", "Rostral_likelihood", "CaudalX", "CaudalY", "Caudal_likelihood"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
[file,path] = uigetfile("*.csv",'Select a DeepLabCut *.csv file');
Position_Table = readtable(fullfile(path,file), opts);

% Convert to output type
position.sample=[];
position.Rostral_XY=[Position_Table.RostralX  Position_Table.RostralY Position_Table.Rostral_likelihood];
position.Caudal_XY=[Position_Table.CaudalX  Position_Table.CaudalY Position_Table.Caudal_likelihood];

% Clear temporary variables
clear opts path;

%% Info about video file (number of frames)

% Import the data
[file,path] = uigetfile("*.avi",'Select a video file for recording');
info = mmfileinfo(fullfile(path,file));
obj = VideoReader(fullfile(path,file));
%disp(['Number of frames: ' num2str(obj.NumFrames), ' - FrameRate: ' num2str(obj.FrameRate)]);

%% Detect rising edge of the TTL AI0 from Intan
read_Intan_RHD2000_file
cd(path)
sample_rate=frequency_parameters.amplifier_sample_rate;
num_channels = length(board_adc_channels); % ADC input info from header file
fileinfo = dir('analogin.dat');
num_samples = fileinfo.bytes/(num_channels * 2); % uint16 = 2 bytes
fid = fopen('analogin.dat', 'r');
v = fread(fid, [num_channels, num_samples], 'uint16');
fclose(fid);
TTL_frame_video = (v-32768) * 0.0003125; % convert to volts

%% sample-frame = alignment
% onset and offset of the TTL
sdnum=0;
meanchannel=mean(TTL_frame_video); 
sdchannel=std(double(TTL_frame_video ));
crossed=TTL_frame_video>meanchannel+(sdnum*sdchannel);
ShiftedData=int16((crossed-circshift(crossed,1)));
TTL_ON_OFF(:,1)=find(ShiftedData==1);
TTL_ON_OFF(:,2)=find(ShiftedData==-1);

clear sdnum meanchannel sdchannel crossed ShiftedData
disp(['Number of TTLs for frames on analog input: ' num2str(length(TTL_ON_OFF))]);
disp(['Number of frames for the video files: ' num2str(obj.NumFrames), ' - FrameRate: ' num2str(obj.FrameRate)]);
disp(['frame difference: #TTLs - #videoFrames: ' num2str(length(TTL_ON_OFF)-obj.NumFrames)])

%if this number is very big there is a problem. Usually there are more TTls
%on analog input as frames in video file. Propably the video recorder stop on PC and send a stop
%signal to camera, but this takes time to arrive and stop the camera FPGA to generate TTLs.

%% create positions
differ_frame_num=(length(TTL_ON_OFF)-obj.NumFrames);
position.sample=TTL_ON_OFF(1:end-differ_frame_num,:);
position.amplifier_sample_rate=sample_rate;
position.video_sample_rate=obj.FrameRate;

for num_frame=1:length(position.sample)
Y1 =  position.Caudal_XY(num_frame,2);
X1=  position.Caudal_XY(num_frame,1);
Y2 = position.Rostral_XY(num_frame,2);
X2 = position.Rostral_XY(num_frame,1);
theta = atan2(X2-X1,Y2-Y1);
position.angle(num_frame,1)=rad2deg(theta);
position.angle(num_frame,2)=(theta);
end

%% Creat transformation matrix -> pixel to mm for the coordinate of the animal
% Constants - size of cage measuerd
monitor_H=500; %cage size in mm 
minitor_W=500;
disp(['cage size is: ' num2str(minitor_W) ' x ' num2str(monitor_H) 'mm - CONSTANT IN THE CODE!!!'])

%% Select corners of the screen
total_frame= round(obj.Duration.*position.video_sample_rate);
obj.CurrentTime=(total_frame./2)./position.video_sample_rate;

disp(['1.) Open video file: ' char(file)])
refImage = readFrame(obj);
disp('2.) Please select the corners of cage (base)!')


fig1=figure(1);
image(refImage);
axis image
hold on
txt = '\leftarrow 1^s^t point';
text(171,47,txt,'FontSize',14,'Color',[1, 0 ,0])

txt = '\leftarrow 2^n^d point';
text(569,42,txt,'FontSize',14,'Color',[1, 0 ,0])

txt = '\leftarrow 3^r^d point';
text(579,443,txt,'FontSize',14,'Color',[1, 0 ,0])

txt = '\leftarrow 4^t^h point';
text(184,447,txt,'FontSize',14,'Color',[1, 0 ,0])

title('draw poly line please follow the order of the cornenr point')

h = drawpolyline(gca,'Color','green');


%% Transformation

movingPoints=h.Position;

fixedPoints=[0,0;minitor_W,0;minitor_W,monitor_H;0,monitor_H];
fixedPoints(:,1)=fixedPoints(:,1)+100; %moving the image in x by 100mm
fixedPoints(:,2)=fixedPoints(:,2)+50; %moving the image in x by 50cm

disp('3.) transformation matrix is DONE!')
tform = fitgeotrans(movingPoints,fixedPoints,'projective'); %calculate transformation matrix

disp('4.) transfrom the image')
Jregistered = imwarp(refImage,tform,'OutputView',imref2d(size(refImage)));


%% Plotting for inspection
disp('5.) plotting figures');
figure;
image(refImage);
xlabel('distance [pixel]')
ylabel('distance [pixel]')
title('original');
axis image;
figure;
image(Jregistered);
xlabel('distance [mm]')
ylabel('distance [mm]')
axis image;
title('cropped and transform image');

close(fig1)
pause(2)
close all

%% Transform pixel to mm

[X_transformd,Y_transformed] = transformPointsForward(tform,Position_Table.RostralX,Position_Table.RostralY);
Position_Table.RostralX_mm = X_transformd;
Position_Table.RostralY_mm = Y_transformed;
position.RostralXY_mm = [X_transformd Y_transformed];


[X_transformd,Y_transformed] = transformPointsForward(tform,Position_Table.CaudalX,Position_Table.CaudalY);
Position_Table.CaudalX_mm = X_transformd;
Position_Table.CaudalY_mm = Y_transformed;
Position_Table.samples_TTL_ON_OFF=position.sample;
position.CaudalXY_mm=[ X_transformd  Y_transformed];
Position_Table.angle_deg_rad=position.angle;

position.tform=tform;
position.filename_video=file;

%% saving variables
save(fullfile(path,file(1:end-4)) ,'position','Position_Table')

%% clear variables
clear all

