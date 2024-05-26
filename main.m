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

sidelobe_tol = 50; % percent of max amplitude

% Simulation config
only_focus_opt = false; % Optimize only focal spots or entire grid
use_greens_fctn = true; % Use Green's function to obtain propagation matrix A (assuming point sources and a lossless homogeneous medium)

get_current_A = "A_2D_3Trs"; % Use precomputed propagation matrix - can be logical or a string containing the file name in Lin_Prop_Matrices
do_time_reversal = false; % Phase retrieval with time reversal as a comparison
save_results = false;

if isempty(dx)
    dx_factor = 1;
else
    dx_factor = (1500 / f0 / 3) / dx; % = (c0 / f0 / ppw) / dx
end

[kgrid, medium, sensor, sensor_mask, b_des, b_des_pl, b_mask, t_mask_ps, karray_t, only_focus_opt, space_limits, ...
    active_ids, mask2el, el_per_t, t_pos, t_rot, ~, ~, point_pos_m, ~, dx_factor1, input_args] = ...
    init(f0, n_dim, dx_factor, ...
    'sidelobe_tol', sidelobe_tol, 't1_scan', t1w_filename, 'ct_scan', ct_filename, 'only_focus_opt', only_focus_opt, 'use_greens_fctn', use_greens_fctn);

[kgridP, mediumP, sensorP, sensor_maskP, ~, ~, ~, t_mask_psP, karray_tP, ~, ~, ...
    active_idsP, ~, ~, ~, ~, plot_offset, point_pos, ~, grid_size, dx_factorP, input_argsP] = ...
    init(f0, n_dim, dx_factor * plot_dx_factor, 't1_scan', t1w_filename, 'ct_scan', ct_filename, 'only_focus_opt', only_focus_opt, 'use_greens_fctn', use_greens_fctn);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Time Reversal
if do_time_reversal
    tr.p = sim_exe(kgrid, medium, sensor, f0, b_des, b_mask, t_mask_ps, false, input_args);
    % tr.p = max(abs(tr.p)) * exp(-1j * angle(tr.p)); % All elements with same amplitude
    tr.p = conj(tr.p); % Var amplitude
    tr.b = sim_exe(kgridP, mediumP, sensorP, f0, tr.p, t_mask_psP, sensor_maskP, true, input_argsP, 'karray_t', karray_tP);
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

obs_ids = reshape(logical(b_mask), numel(b_mask), 1);
opt_ids = limit_space(medium.sound_speed);

if only_focus_opt
    % Preserve sonicated points only
    ip.A = A(obs_ids, :);
    b_ip_des = b_des; % = b_des_pl(obs_ids)

    init_ids = true(size(ip.A, 1), 1);
    beta_L2 = 0.0;
else
    % Take entire observation grid into account
    ip.A = A;
    b_ip_des = b_des_pl;

    init_ids = get_init_ids(kgrid, min(medium.sound_speed(:)) / f0, b_mask);
    beta_L2 = 0.0;
end

p_init = pinv(ip.A(init_ids, :)) * b_ip_des(init_ids, :);

ip.p = solvePhasesOnly(ip.A, b_ip_des, opt_ids, p_init, init_ids, beta_L2, mask2el, el_per_t, true); % Amp fixed
ip.p_gt = solvePhases_Amp(ip.A, b_ip_des, opt_ids, p_init, init_ids, beta_L2); % var Amp

ip.t_solve = toc;

%% Obtain Acoustic profile
% ip.b_gt = sim_exe(kgridP, mediumP, sensorP, f0, ip.p, t_mask_psP, sensor_maskP, true, input_argsP, 'karray_t', karray_tP);
ip.b = A * ip.p;
ip.b = reshape(ip.b, size(kgrid.k));

ip.b_gt = A * ip.p_gt;
ip.b_gt = reshape(ip.b_gt, size(kgrid.k));

%% Save Results in mat-file
current_datetime = string(datestr(now, 'yyyymmddHHMMSS'));
if save_results
    res_filename = "results";
    if ~only_focus_opt
        ip.A = []; % A might be very large...
    end
    save(fullfile("Results", current_datetime + "_" + res_filename + ".mat"), ...
        "f0", "kgrid", "b_mask", "t_mask_ps", "active_ids", "init_ids", "mask2el", "t_pos", "t_rot", "tr", "ip", "b_ip_des", "point_pos", "point_pos_m", ...
        "only_focus_opt", "input_args");
end

%% TR Results
if do_time_reversal
    plot_results(kgridP, tr.p, tr.b, 'Time Reversal', mask2el, t1w_filename, plot_offset, grid_size, dx_factorP, save_results, current_datetime, 'slice', point_pos.slice);
end

%% IP Results
plot_results(kgrid, ip.p, ip.b, 'Inverse Problem', mask2el, t1w_filename, plot_offset, grid_size, dx_factor1, save_results, current_datetime, 'slice', point_pos.slice);
plot_results(kgrid, ip.p_gt, ip.b_gt, 'IP2', mask2el, t1w_filename, plot_offset, grid_size, dx_factor1, save_results, current_datetime, 'slice', point_pos.slice);
% plot_results(kgridP, ip.p, ip.b_gt, 'Ground Truth', mask2el, t1w_filename, plot_offset, grid_size, dx_factor1, save_results, current_datetime, 'slice', point_pos.slice);

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

if only_focus_opt
    err_ip = real_ip - b_ip_des;
    err_ip_gt = real_ip_gt - b_ip_des;
    
    std_ip = sqrt(sum(err_ip.^2) / length(err_ip));
    std_ip_gt = sqrt(sum(err_ip_gt.^2) / length(err_ip_gt));
    
    maxDev_ip = max(abs(err_ip));
    maxDev_ip_gt = max(abs(err_ip_gt));
    
    % fprintf("\nDesired Amplitudes (kPa):\n")
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
else
    offTar_real_ip = real_ip;
    offTar_real_ip(obs_ids | ~opt_ids) = [];
    offTar_real_ip_gt = real_ip_gt;
    offTar_real_ip_gt(obs_ids | ~opt_ids) = [];

    tar_real_ip = real_ip(obs_ids);
    tar_real_ip_gt = real_ip_gt(obs_ids);

    init_real_ip = real_ip(init_ids);
    init_real_ip_gt = real_ip_gt(init_ids);

    fprintf("\nInverse Problem Init Points (kPa):\n")
    disp(init_real_ip' * 1e-3)
    disp(init_real_ip_gt' * 1e-3)

    fprintf("\nInverse Problem Min Target Area (kPa):\n")
    disp(min(tar_real_ip) * 1e-3)
    disp(min(tar_real_ip_gt) * 1e-3)

    fprintf("\nInverse Problem Max Off-Target (kPa):\n")
    disp(max(offTar_real_ip) * 1e-3)
    disp(max(offTar_real_ip_gt) * 1e-3)

    point_coord = [-16, -27];
    point_val_ip = abs(ip.b(plot_offset(1) + point_coord(1), plot_offset(3) + point_coord(2)));
    point_val_ip_gt = abs(ip.b_gt(plot_offset(1) + point_coord(1), plot_offset(3) + point_coord(2)));

    fprintf("\nInverse Problem sel. Point Pressure (kPa):\n")
    disp(max(point_val_ip) * 1e-3)
    disp(max(point_val_ip_gt) * 1e-3)
end
