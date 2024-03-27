clear;

%% Load desired Results file
in_filename = "20240325190245_results";
load(fullfile("..", "Results", in_filename + ".mat"));

%% Compute Delays for Time Reversal and Inverse Problem
tr_angles = angle(tr.p') + pi; % between 0 and 2*pi
tr_delays = tr_angles / (2*pi * f0);

ip_angles = angle(ip.p') + pi; % between 0 and 2*pi
ip_delays = ip_angles / (2*pi * f0);

%% Save as txt files
in_filename = char(in_filename);
record_datetime = string(in_filename(1:14));
out_filename = "delayParam";

dlmwrite(fullfile("..", "Transducer_Delay_Files", record_datetime + "_" + out_filename + "_tr.txt"), tr_delays, 'delimiter', ' ');
dlmwrite(fullfile("..", "Transducer_Delay_Files", record_datetime + "_" + out_filename + "_ip.txt"), ip_delays, 'delimiter', ' ');

