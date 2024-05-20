clear;
close all;

%% Init
f0 = 470e3; % Hz - transducer frequency
n_dim = 2;
dx = 1e-3;
% dx = [];
plot_dx_factor = 1;

t1w_filename = fullfile('Scans', 'dummy_t1w.nii');
ct_filename = fullfile('Scans', 'dummy_pseudoCT.nii');
t1w_filename = [];
ct_filename = [];

% Simulation config
only_focus_opt = false; % Optimize only focal spots or entire grid
use_greens_fctn = true; % Use Green's function to obtain propagation matrix A (assuming point sources and a lossless homogeneous medium)

get_current_A = "A_2D_3Trs"; % Use precomputed propagation matrix - can be logical or a string containing the file name in Lin_Prop_Matrices
get_current_AP = false; % Use precomputed propagation matrix - Only to plot resulting acoustic profile
do_time_reversal = false; % Phase retrieval with time reversal as a comparison
save_results = false;

if isempty(dx)
    dx_factor = 1;
else
    dx_factor = (1500 / f0 / 3) / dx; % = (c0 / f0 / ppw) / dx
end

[kgrid, medium, sensor, sensor_mask, b_des, b_des_pl, b_mask, t_mask_ps, karray_t, only_focus_opt, space_limits, ...
    active_ids, mask2el, el_per_t, t_pos, t_rot, amp_in, ~, ~, point_pos_m, ~, ~, input_args] = ...
    init(f0, n_dim, dx_factor, 't1_scan', t1w_filename, 'ct_scan', ct_filename, 'only_focus_opt', only_focus_opt, 'use_greens_fctn', use_greens_fctn);

[kgridP, mediumP, sensorP, sensor_maskP, ~, ~, ~, t_mask_psP, karray_tP, ~, ~, ...
    active_idsP, ~, ~, ~, ~, ~, plot_offset, point_pos, ~, grid_size, dx_factorP, input_argsP] = ...
    init(f0, n_dim, dx_factor * plot_dx_factor, 't1_scan', t1w_filename, 'ct_scan', ct_filename, 'only_focus_opt', only_focus_opt, 'use_greens_fctn', use_greens_fctn);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Time Reversal
if do_time_reversal
    tr.p = sim_exe(kgrid, medium, sensor, f0, b_des, b_mask, t_mask_ps, false, input_args);
    tr.p = max(abs(tr.p)) * exp(-1j * angle(tr.p)); % All elements with same amplitude
    tr.b = sim_exe(kgrid, medium, sensor, f0, tr.p, t_mask_ps, sensor_mask, true, input_args, 'karray_t', karray_t);
else
    tr = [];
end

%% Obtain propagation operator -> acousticFieldPropagator (Green's functions)
% A = linearPropagator_vs_acousticFieldPropagator(t_mask, f0, medium.sound_speed, kgrid.dx);
if ~exist('A', 'var')
    A = obtain_linear_propagator(kgrid, medium, sensor, sensor_mask, input_args, t_mask_ps, karray_t, f0, get_current_A, use_greens_fctn, ...
        'active_ids', active_ids);
end

%% Solve Inverse Problem
tic
obs_ids = find(b_mask);
if true%only_focus_opt
    % Preserve sonicated points only
    ip.A = A(obs_ids, :);
    b_ip_des = b_des;

    activeA_ids = 1:size(ip.A, 1);

%     ip.sq_beta = 0;%0.1; % constraint scaling factor
else
    % Take entire observation grid into account
    ip.A = A;
    b_ip_des = b_des_pl;

    activeA_ids = zeros(1, size(A, 1));
    activeA_ids(obs_ids) = 1;
    activeA_ids = logical(activeA_ids);

%     [ip.A, b_ip_des, activeA_ids] = limit_space(ip.A, b_ip_des, activeA_ids, space_limits, plot_offset, dx_factor, medium.sound_speed);

%     ip.sq_beta = 0;%100; % constraint scaling factor
end

% p_init = pinv(ip.A(activeA_ids, :)) * b_ip_des(activeA_ids, :);
beta_L2 = 0.1;

ip.p_gt = solvePhases_Amp(ip.A, b_ip_des, p_init, beta_L2); % var Amp

ip.p = solvePhasesOnly(ip.A(activeA_ids, :), ip.A(~activeA_ids, :), b_ip_des(activeA_ids, :), max(b_ip_des) / 10, p_init, beta_L2, mask2el, el_per_t, true); % Amp fixed

% ip.p_gt = solvePhasesOnly(ip.A(activeA_ids, :), ip.A(~activeA_ids, :), b_ip_des(activeA_ids, :), max(b_ip_des) / 10, mask2el, el_per_t, false); % Amp fixed

ip.t_solve = toc;


% Evaluate obtained phase terms in forward simulation
if plot_dx_factor ~= 1
    clear A;
    A = obtain_linear_propagator(kgridP, mediumP, sensorP, sensor_maskP, input_argsP, t_mask_psP, karray_tP, f0, get_current_AP, use_greens_fctn, ...
        'active_ids', active_idsP); % Obtain high resolution A - Discard if no point sources!
end

% ip.b_gt = sim_exe(kgridP, mediumP, sensorP, f0, ip.p, t_mask_psP, sensor_maskP, true, input_argsP, 'karray_t', karray_tP);
ip.b = A * ip.p;
ip.b = reshape(ip.b, size(kgridP.k));

ip.b_gt = A * ip.p_gt;
ip.b_gt = reshape(ip.b_gt, size(kgridP.k));

%% Save Results in mat-file
current_datetime = string(datestr(now, 'yyyymmddHHMMSS'));
if save_results
    res_filename = "results";
    if ~only_focus_opt
        ip.A = []; % A might be very large...
    end
    save(fullfile("Results", current_datetime + "_" + res_filename + ".mat"), ...
        "f0", "kgrid", "b_mask", "t_mask_ps", "active_ids", "mask2el", "t_pos", "t_rot", "tr", "ip", "amp_in", "point_pos", "point_pos_m", ...
        "only_focus_opt", "input_args");
end

%% TR Results
if do_time_reversal
    plot_results(kgridP, tr.p, tr.b, 'Time Reversal', mask2el, t1w_filename, plot_offset, grid_size, dx_factorP, save_results, current_datetime);
end

%% IP Results
plot_results(kgridP, ip.p, ip.b, 'Inverse Problem', mask2el, t1w_filename, plot_offset, grid_size, dx_factorP, save_results, current_datetime, 'slice', point_pos.slice);
plot_results(kgridP, ip.p_gt, ip.b_gt, 'Inverse Problem 2', mask2el, t1w_filename, plot_offset, grid_size, dx_factorP, save_results, current_datetime, 'slice', point_pos.slice);
% plot_results(kgridP, ip.p, ip.b_gt, 'Ground Truth', mask2el_delayFiles, t1w_filename, plot_offset, grid_size, dx_factorP, save_results, current_datetime, 'slice', point_pos.slice);

% err = abs(ip.b) - abs(ip.b_gt);
% plot_results(kgridP, [], err, 'Difference', mask2el, t1w_filename, plot_offset, grid_size, dx_factorP, save_results, current_datetime, 'slice', point_pos.slice);
% figure
% histogram(err(:))
% xlabel("Pressure Deviation (Pa)")

%% Metrics evaluation
disp("Time until solver converged: " + string(ip.t_solve) + " s")
% fprintf("\nDesired Transducer Amplitude (kPa):\n")
% disp(ip.u(1) * 1e-3)
% fprintf("\nTransducer Amplitude mean deviation (Pa):\n")
% disp(mean(abs(abs(ip.p) - ip.u)))


real_ip = abs(ip.A * ip.p);
real_ip_gt = abs(ip.A * ip.p_gt);

err_ip = real_ip - b_ip_des;
err_ip_gt = real_ip_gt - b_ip_des;

std_ip = std(err_ip);
std_ip_gt = std(err_ip_gt);

maxDev_ip = max(abs(err_ip));
maxDev_ip_gt = max(abs(err_ip_gt));

% fprintf("\nInput Amplitudes (kPa):\n")
% disp(b_ip_des' * 1e-3)
% 
% fprintf("\nReal Amplitudes (kPa):\n")
% disp(real_ip' * 1e-3)
% disp(real_ip_gt' * 1e-3)

fprintf("\nInverse Problem STD (kPa):\n")
disp(std_ip * 1e-3)
disp(std_ip_gt * 1e-3)

fprintf("\nInverse Problem Max Dev (kPa):\n")
disp(maxDev_ip * 1e-3)
disp(maxDev_ip_gt * 1e-3)

