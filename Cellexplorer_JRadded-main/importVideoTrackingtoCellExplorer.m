sr = 20000.0; 
criteria = 0.95;
gaussian_kernel = 20;
basename = 'amplifier';

rostral_x_mm = position.RostralXY_mm(:,1);
rostral_y_mm = position.RostralXY_mm(:,2);
rostral_likelihood = position.Rostral_XY(:,3);

caudal_x_mm = position.CaudalXY_mm(:,1);
caudal_y_mm = position.CaudalXY_mm(:,2);
caudal_likelihood = position.Caudal_XY(:,3);

time = position.sample(:,1)/sr;

good_rostral = rostral_likelihood > criteria;
good_caudal = caudal_likelihood > criteria; 
good_position = and(good_rostral, good_caudal);

good_rostral_x_mm = rostral_x_mm(good_position);
good_rostral_y_mm = rostral_y_mm(good_position);
good_caudal_x_mm = caudal_x_mm(good_position);
good_caudal_y_mm = caudal_y_mm(good_position);
good_time = time(good_position);
good_angle = position.angle(good_position,1);

smooth_good_rostral_x_mm = smoothdata(good_rostral_x_mm,'gaussian',gaussian_kernel);
smooth_good_rostral_y_mm = smoothdata(good_rostral_y_mm,'gaussian',gaussian_kernel);

dt = diff(good_time);
dx = diff(smooth_good_rostral_x_mm);
dy = diff(smooth_good_rostral_y_mm);
dr = sqrt(dx.^2 + dy.^2);
rostral_speed = dr ./ dt;

video_tracking = struct;
video_tracking.timestamps = good_time;
video_tracking.sr = position.video_sample_rate;

video_tracking.position = struct; 
video_tracking.position.x = smooth_good_rostral_x_mm ./ 10; 
video_tracking.position.y = smooth_good_rostral_y_mm ./ 10; 
video_tracking.position.units = 'centimeters';

video_tracking.speed = rostral_speed; 
video_tracking.time_series = good_angle;
video_tracking.notes = 'time_series: head direction angle in degrees, as estimated from the video using the two LEDs';

save(strcat(basename,'.video_tracking.behavior.mat'),'video_tracking');