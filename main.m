clear;
close all;

%% Init
f0 = 500e3; % Hz - transducer frequency
n_dim = 2;
dx_factor = 1;
if n_dim == 2
    grid_size = [192, 256] * 1e-3; % m in [x, y] respectively
else
    grid_size = [120, 120, 100] * 1e-3; % m in [x, y, z] respectively
end
[kgrid, medium, ppp] = init_grid_medium(f0, grid_size, 'n_dim', n_dim, 'dx_factor', 1 / dx_factor);
[sensor, sensor_mask] = init_sensor(kgrid, ppp);

% Scan init
t1w_filename = fullfile('Scans', 'dummy_t1w.nii');
t1w_offset = [96, 127, 126] + 1; % Offset to Scan center

% Simulation config
only_focus_opt = true; % Optimize only focal spots or entire grid
use_greens_fctn = true; % Use Green's function to obtain propagation matrix A (assuming point sources and a lossless homogeneous medium)

get_current_A = false; % Use precomputed propagation matrix - can be logical or a string containing the file name in Lin_Prop_Matrices
do_time_reversal = false; % Phase retrieval with time reversal as a comparison
save_results = false;

%% Define Transducer Geometry

if kgrid.dim == 2
    % Linear array
    t1_pos = [-63, -2]';
    t2_pos = [0, 68]';
    t3_pos = [73, -2]';
    t_pos = [t1_pos, t2_pos, t3_pos];
    t_rot = [false, true, false];

    num_elements = 50;
    spacing = ceil(dx_factor);

    t_mask_ps = zeros(kgrid.Nx, kgrid.Ny);
    for i = 1:length(t_rot)
        el_offset = round((t1w_offset(1) + t_pos(1, i)) * dx_factor); % grid points
        shift = round((t1w_offset(3) + t_pos(2, i)) * dx_factor); % tangential shift in grid points
    
        t_mask_ps = t_mask_ps + create_linear_array(kgrid, num_elements, el_offset, shift, spacing, t_rot(i));
    end
    

    karray_t = [];
    active_ids = [];
    mask2el_delayFiles = [];
    t_mask_ps = t_mask_ps > 0; % Return to logical in case of overlaps

%     imagesc(t_mask_ps, [-1 1])
%     colormap(getColorMap);
else
    % Planar Array
    if use_greens_fctn
        t_name = "std";
    else
        t_name = "std_orig";
    end
    sparsity_name = "sparsity_ids";
    t1_pos = [-55, 20, 0]' * 1e-3; % m
    t1_rot = [-90, 0, 90]'; % deg
    t2_pos = [20, -55, 0]' * 1e-3; % m
    t2_rot = [-90, 0, 180]'; % deg

    t_pos = [t1_pos, t2_pos];
    t_rot = [t1_rot, t2_rot];
    active_tr_ids = [1, 2];

    [karray_t, t_mask_ps, active_ids, mask2el_delayFiles] = create_transducer(kgrid, t_name, sparsity_name, t_pos, t_rot, active_tr_ids);

%     voxelPlot(double(t_mask))
end

%% Define (intracranial) Beamforming Pattern

if kgrid.dim == 2

    % Focal points - rel. to transducer surface
    point_pos_m.x = [-16, 23];
    point_pos.slice = 32;
    point_pos_m.y = [-27, -26];
    amp_in = [200, 200]' * 1e3; % Pa

    point_pos.x = round((t1w_offset(1) + point_pos_m.x) * dx_factor);
    point_pos.y = round((t1w_offset(3) + point_pos_m.y) * dx_factor);

    % Assign amplitude acc. to closest position
    idx = sub2ind([kgrid.Nx, kgrid.Ny], point_pos.x, point_pos.y);
    [~, order] = sort(idx);
    amp_in = amp_in(order);

    b_mask = zeros(kgrid.Nx, kgrid.Ny);
    
    if ~only_focus_opt

        % Stimulate Disc pattern
        for i = 1:length(point_pos.x)
            b_mask = b_mask + makeDisc(kgrid.Nx, kgrid.Ny, point_pos.x(i), point_pos.y(i), round(0.025 * kgrid.Nx), false);
            amp_in = amp_in(i) * ones(sum(b_mask(:)), 1);
        end

    else

        for point = 1:length(point_pos.x)
            b_mask(point_pos.x(point), point_pos.y(point)) = 1;
        end
    end

    f = figure;
    f.Position = [700 485 484 512];
    imagesc(imrotate(b_mask + t_mask_ps, 90), [-1 1])
    colormap(getColorMap);
    title("Setup")
else
    only_focus_opt = true;

    % Focal points - rel. to transducer surface
    point_pos_m.x = [75] * 1e-3; % m
    point_pos_m.y = [0] * 1e-3; % m
    point_pos_m.z = [0] * 1e-3; % m
    amp_in = [100]' * 1e3; % Pa

    point_pos.x = point_pos_m.x + t_pos(1, 1);
    point_pos.y = point_pos_m.y + t_pos(2, 1);
    point_pos.z = point_pos_m.z + t_pos(3, 1);

    point_pos.x = round((point_pos.x - kgrid.x_vec(1)) / kgrid.dx); % grid points
    point_pos.y = round((point_pos.y - kgrid.y_vec(1)) / kgrid.dy); % grid points
    point_pos.z = round((point_pos.z - kgrid.z_vec(1)) / kgrid.dz); % grid points

    % Assign amplitude acc. to closest position
    idx = sub2ind([kgrid.Nx, kgrid.Ny, kgrid.Nz], point_pos.x, point_pos.y, point_pos.z);
    [~, order] = sort(idx);
    amp_in = amp_in(order);

    b_mask = zeros(kgrid.Nx, kgrid.Ny, kgrid.Nz);

    for point = 1:length(point_pos.x)
        b_mask(point_pos.x(point), point_pos.y(point), point_pos.z(point)) = 1;
    end

    voxelPlot(double(t_mask_ps + b_mask))
end

% Create desired signal
phase = zeros(length(amp_in), 1); % Zero phase for entire observation plane

b_des = amp_in .* exp(1j*phase); % only observed elements

b_des_pl = zeros(kgrid.Nx * kgrid.Ny, 1); % entire plane
b_des_pl(find(b_mask)) = b_des;

% set simulation input options
input_args = {'PMLSize', 'auto', 'PMLInside', false, 'PlotPML', true, 'DisplayMask', b_mask + t_mask_ps, 'RecordMovie', false};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Time Reversal
if do_time_reversal
    tr.p = sim_exe(kgrid, medium, sensor, f0, b_des, b_mask, t_mask_ps, false, input_args);
    tr.p = max(abs(tr.p)) * exp(-1j * angle(tr.p)); % All elements with same amplitude
    tr.b = sim_exe(kgrid, medium, sensor, f0, tr.p, t_mask_ps, sensor_mask, true, input_args, 'karray_t', karray_t);
else
    tr = [];
end

%% Inverse Problem

% Obtain propagation operator -> acousticFieldPropagator (Green's functions)
% A = linearPropagator_vs_acousticFieldPropagator(t_mask, f0, medium.sound_speed, kgrid.dx);
A = obtain_linear_propagator(kgrid, medium, sensor, sensor_mask, input_args, t_mask_ps, karray_t, f0, get_current_A, use_greens_fctn, ...
    'active_ids', active_ids);

% Solve inverse problem
tic
obs_ids = find(b_mask);
if only_focus_opt
    % Preserve sonicated points only
    ip.A = A(obs_ids, :);
    b_ip_des = b_des;

%     ip.sq_beta = 0;%0.1; % constraint scaling factor
else
    % Take entire observation grid into account
    ip.A = A;
    b_ip_des = b_des_pl;

%     ip.sq_beta = 0;%100; % constraint scaling factor
end

ip.p = pinv(ip.A) * b_ip_des; % Initial solution cosidering phases

% ip.u = max(abs(ip.p)) * ones(size(ip.p)); % universal transducer amplitude
% 
% % Apply constraints (same transducer amplitude among elements)
% ip.A = [ip.A; ip.sq_beta * eye(length(ip.p))];
% b_ip_des = [b_ip_des; ip.sq_beta * ip.u];

% Solve phase retrieval problem
opts = struct;
opts.initMethod = 'custom';
opts.customx0 = ip.p;

[ip.p, outs, opts] = solvePhaseRetrieval(ip.A, ip.A', b_ip_des, [], opts); % var Amplitude

ip.p = solvePhasesOnly(A(obs_ids, :), A(~obs_ids, :), b_des, max(b_des) / 10); % Amplitude fixed

ip.t_solve = toc;

% ip.p = max(abs(ip.p)) * exp(1j * angle(ip.p)); % All elements with same amplitude

% Evaluate obtained phase terms in forward simulation

% ip.b_gt = sim_exe(kgrid, medium, sensor, f0, ip.p, t_mask_ps, sensor_mask, true, input_args, 'karray_t', karray_t);
ip.b = A * ip.p;
ip.b = reshape(ip.b, size(kgrid.k));

%% Save Results in mat-file
if save_results
    current_datetime = string(datestr(now, 'yyyymmddHHMMSS'));
    res_filename = "results";
    if ~only_focus_opt
        ip.A = []; % A might be very large...
    end
    save(fullfile("Results", current_datetime + "_" + res_filename + ".mat"), ...
        "f0", "kgrid", "b_mask", "t_mask_ps", "active_ids", "mask2el_delayFiles", "t_pos", "t_rot", "tr", "ip", "amp_in", "point_pos", "point_pos_m", ...
        "only_focus_opt", "input_args");
end

%% TR Results
if do_time_reversal
    plot_results(kgrid, tr.p, tr.b, t_pos, 'Time Reversal', t1w_filename, t1w_offset);
end

%% IP Results
plot_results(kgrid, ip.p, ip.b, t_pos, 'Inverse Problem', t1w_filename, t1w_offset, 'slice', point_pos.slice);

% Metrics evaluation
disp("Time until solver converged: " + string(ip.t_solve) + " s")
fprintf("\nDesired Transducer Amplitude (kPa):\n")
disp(ip.u(1) * 1e-3)
fprintf("\nTransducer Amplitude mean deviation (Pa):\n")
disp(mean(abs(abs(ip.p) - ip.u)))

% Evaluation of each defined point
if only_focus_opt
    b_tr_points = [];
    b_ip_points = [];

    if kgrid.dim == 2
        for point = 1:length(point_pos.x)
%             b_tr_points = [b_tr_points, tr.b(point_pos.x(point), point_pos.y(point))];
            b_ip_points = [b_ip_points, ip.b(point_pos.x(point), point_pos.y(point))];
        end
    else
        for point = 1:length(point_pos.x)
%             b_tr_points = [b_tr_points, tr.b(point_pos.x(point), point_pos.y(point), point_pos.z(point))];
            b_ip_points = [b_ip_points, ip.b(point_pos.x(point), point_pos.y(point), point_pos.z(point))];
        end
    end
    
    fprintf("\nInput Amplitudes (kPa):\n")
    disp(amp_in' * 1e-3)
%     fprintf("\nTime Reversal Total Amplitudes (kPa):\n")
%     disp(abs(b_tr_points) * 1e-3)
    fprintf("\nInverse Problem Total Amplitudes (kPa):\n")
    disp(abs(b_ip_points) * 1e-3)
%     fprintf("\nTime Reversal Phase Angles (deg):\n")
%     disp(angle(b_tr_points) / pi * 180)
%     fprintf("\nInverse Problem Phase Angles (deg):\n")
%     disp(angle(b_ip_points) / pi * 180)
end

