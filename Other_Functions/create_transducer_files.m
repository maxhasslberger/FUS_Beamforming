clear;

%% Load desired Results file
in_filename = "20240327204947_results";
load(fullfile("..", "Results", in_filename + ".mat"));

for i = 1:size(mask2el_ids, 2) % for each transducer

    %% Compute Delays for Inverse Problem
    ip_corrected = ip.p(mask2el_ids(:, i));
    ip_angles = angle(ip_corrected') + pi; % between 0 and 2*pi
    ip_delays = ip_angles / (2*pi * f0);
    
    %% Save as txt files
    in_filename = char(in_filename);
    record_datetime = string(in_filename(1:14));
    out_filename = "tr" + string(i);
    
    dlmwrite(fullfile("..", "Transducer_Delay_Files", record_datetime + "_" + out_filename + ".txt"), ip_delays, 'delimiter', ' ');
end
