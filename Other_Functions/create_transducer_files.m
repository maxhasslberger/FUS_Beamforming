clear;

%% Load desired Results file
in_filename = "20240325164602_results";
load(fullfile("..", "Results", in_filename + ".mat"));

%% Compute Delays for Time Reversal and Inverse Problem
tr_angles = angle(tr.p') + pi; % between 0 and 2*pi
tr_delays = tr_angles / (2*pi * f0);

ip_angles = angle(ip.p') + pi; % between 0 and 2*pi
ip_delays = ip_angles / (2*pi * f0);

%% Save as txt files
current_datetime = string(datestr(now, 'yyyymmddHHMMSS'));
out_filename = "test";

dlmwrite(fullfile("..", "Transducer_Delay_Files", current_datetime + "_" + out_filename + "_tr.txt"), tr_delays, 'delimiter', ' ');
dlmwrite(fullfile("..", "Transducer_Delay_Files", current_datetime + "_" + out_filename + "_ip.txt"), ip_delays, 'delimiter', ' ');

