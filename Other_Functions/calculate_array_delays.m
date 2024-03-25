clear;

%% Load desired Results file
in_filename = "res";
load(fullfile("..", in_filename + ".mat"));

%% Compute Delays for Time Reversal and Inverse Problem
tr_angles = angle(tr.p_tr') + pi; % between 0 and 2*pi
tr_delays = tr_angles / (2*pi * f0);

ip_angles = angle(ip.p_ip') + pi; % between 0 and 2*pi
ip_delays = ip_angles / (2*pi * f0);

%% Save as txt files
current_datetime = string(datestr(now, 'yyyymmddHHMMSS'));

dlmwrite(fullfile("..", "Transducer_Delay_Files", current_datetime + "_tr.txt"), tr_delays, 'delimiter', ' ');
dlmwrite(fullfile("..", "Transducer_Delay_Files", current_datetime + "_ip.txt"), ip_delays, 'delimiter', ' ');

