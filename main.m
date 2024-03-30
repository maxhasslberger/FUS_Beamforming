clear;
close all;

%% Init

f0 = 500e3; % Hz - transducer frequency
n_dim = 2;
dx_factor = 1;
if n_dim == 2
    grid_size = [100, 100] * 1e-3; % m in [x, y] respectively
else
    grid_size = [120, 120, 100] * 1e-3; % m in [x, y, z] respectively
end
[kgrid, medium, ppp] = init_grid_medium(f0, grid_size, 'n_dim', n_dim, 'dx_factor', 1 / dx_factor);
[sensor, sensor_mask] = init_sensor(kgrid, ppp);

only_focus_opt = true; % Optimize only focal spots or entire grid
set_current_A = false; % Use precomputed propagation matrix - can be logical or a string containing the file name in Lin_Prop_Matrices
do_time_reversal = true;
save_results = true;

%% Define Transducer Geometry

if kgrid.dim == 2
    % Linear array
    t1_pos = [-45, 20]' * 1e-3; % m
    t2_pos = [-20, -45]' * 1e-3; % m
    t_pos = [t1_pos, t2_pos];
    t_rot = [];

    el1_offset = round((t1_pos(1) - kgrid.x_vec(1)) / kgrid.dx); % grid points
    el2_offset = round((t2_pos(2) - kgrid.y_vec(1)) / kgrid.dy); % grid points
    shift1 = round(t1_pos(2) / kgrid.dx); % m -> tangential shift in grid points
    shift2 = round(t2_pos(1) / kgrid.dx); % m -> tangential shift in grid points
    if only_focus_opt
        num_elements = 50;
        spacing = ceil(1e-3 / kgrid.dx * dx_factor); % m -> grid points between elements

        t_mask = create_linear_array(kgrid, num_elements, el1_offset, shift1, spacing, false);
        % t_mask = t_mask + create_linear_array(kgrid, num_elements, el2_offset, shift2, spacing, true); % Second (orthogonal) linear array
        % t_mask = t_mask + create_linear_array(kgrid, num_elements, round(kgrid.Nx - el1_offset), shift1, spacing, false); % Second (antiparallel) linear array
        % t_mask = t_mask + create_linear_array(kgrid, num_elements, round(kgrid.Nx - el2_offset), shift2, spacing, true); % Third (orthogonal) linear array
    else
        num_elements = 25;
        spacing = ceil(2e-3 / kgrid.dx * dx_factor); % m -> grid points between elements

        t_mask = create_linear_array(kgrid, num_elements, el1_offset, shift1, spacing, false);
        t_mask = t_mask + create_linear_array(kgrid, num_elements, el2_offset, shift2, spacing, true); % Second (orthogonal) linear array
        % t_mask = t_mask + create_linear_array(kgrid, num_elements, round(kgrid.Nx - el1_offset), shift1, spacing, false); % Second (antiparallel) linear array
        % t_mask = t_mask + create_linear_array(kgrid, num_elements, round(kgrid.Nx - el2_offset), shift2, spacing, true); % Third (orthogonal) linear array
    end
    
    karray_t = [];
    mask_active = [];
    mask2el_delayFiles = [];
    t_mask = t_mask > 0; % Return to logical in case of overlaps

%     imagesc(t_mask, [-1 1])
%     colormap(getColorMap);
else
    % Planar Array
    t_name = "std";
    sparsity_name = "sparsity_ids";
    t1_pos = [-55, 20, 0]' * 1e-3; % m
    t1_rot = [-90, 0, 90]'; % deg
    t2_pos = [20, -55, 0]' * 1e-3; % m
    t2_rot = [-90, 0, 180]'; % deg

    t_pos = [t1_pos, t2_pos];
    t_rot = [t1_rot, t2_rot];
    active_tr_ids = [1];

    [karray_t, mask_active, mask2el_delayFiles] = create_transducer(kgrid, t_name, sparsity_name, t_pos, t_rot, active_tr_ids);
    t_mask = karray_t.getArrayBinaryMask(kgrid);

%     voxelPlot(double(t_mask))
end

%% Define (intracranial) Beamforming Pattern

if kgrid.dim == 2
    
    if ~only_focus_opt

        % Ring
        b_mask = makeDisc(kgrid.Nx, kgrid.Ny, round(0.7 * kgrid.Nx), round(0.7 * kgrid.Ny), round(0.2 * kgrid.Nx), false) ...
            - makeDisc(kgrid.Nx, kgrid.Ny, round(0.7 * kgrid.Nx), round(0.7 * kgrid.Ny), round(0.15 * kgrid.Nx), false);
        amp_in = 30000 * ones(sum(b_mask(:)), 1);

        point_pos_m = [];
        point_pos = [];
    else

        % Points
        point_pos_m.x = [10, -10] * 1e-3; % m
        point_pos_m.y = [0, 35] * 1e-3; % m
        amp_in = [200, 200]' * 1e3; % Pa

        point_pos.x = round((point_pos_m.x - kgrid.x_vec(1)) / kgrid.dx); % grid points
        point_pos.y = round((point_pos_m.y - kgrid.y_vec(1)) / kgrid.dy); % grid points
    
        % Assign amplitude acc. to closest position
        idx = sub2ind([kgrid.Nx, kgrid.Ny], point_pos.x, point_pos.y);
        [~, order] = sort(idx);
        amp_in = amp_in(order);
    
        b_mask = zeros(kgrid.Nx, kgrid.Ny);

        for point = 1:length(point_pos.x)
            b_mask(point_pos.x(point), point_pos.y(point)) = 1;
        end
    end

    imagesc(b_mask + t_mask, [-1 1])
    colormap(getColorMap);
    title("Setup")
else
    only_focus_opt = true;

    % Point - rel. to transducer surface -> [40, 10, 25] mm - TODO: Output point coordinates relative to transducer 1 surface
    point_pos_m.x = [20] * 1e-3; % m
    point_pos_m.y = [0] * 1e-3; % m
    point_pos_m.z = [10] * 1e-3; % m
    amp_in = [100]' * 1e3; % Pa

    point_pos.x = round((point_pos_m.x - kgrid.x_vec(1)) / kgrid.dx); % grid points
    point_pos.y = round((point_pos_m.y - kgrid.y_vec(1)) / kgrid.dy); % grid points
    point_pos.z = round((point_pos_m.z - kgrid.z_vec(1)) / kgrid.dz); % grid points

    % Assign amplitude acc. to closest position
    idx = sub2ind([kgrid.Nx, kgrid.Ny, kgrid.Nz], point_pos.x, point_pos.y, point_pos.z);
    [~, order] = sort(idx);
    amp_in = amp_in(order);

    b_mask = zeros(kgrid.Nx, kgrid.Ny, kgrid.Nz);

    for point = 1:length(point_pos.x)
        b_mask(point_pos.x(point), point_pos.y(point), point_pos.z(point)) = 1;
    end

    voxelPlot(double(t_mask + b_mask))
end

% Create desired signal
phase = zeros(length(amp_in), 1); % Zero phase for entire observation plane

b_des = amp_in .* exp(1j*phase); % only observed elements

b_des_pl = zeros(kgrid.Nx * kgrid.Ny, 1); % entire plane
b_des_pl(find(b_mask)) = b_des;

% set simulation input options
input_args = {'PMLSize', 'auto', 'PMLInside', false, 'PlotPML', true, 'DisplayMask', b_mask + t_mask, 'RecordMovie', false};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Time Reversal
if do_time_reversal
    tr.p = sim_exe(kgrid, medium, sensor, f0, b_des, b_mask, t_mask, false, input_args);
    tr.p = max(abs(tr.p)) * exp(-1j * angle(tr.p)); % All elements with same amplitude
    tr.b = sim_exe(kgrid, medium, sensor, f0, tr.p, t_mask, sensor_mask, true, input_args, 'karray_t', karray_t);
else
    tr = [];
end

%% Inverse Problem

% Obtain propagation operator -> acousticFieldPropagator (Green's functions)
% A = linearPropagator_vs_acousticFieldPropagator(t_mask, f0, medium.sound_speed, kgrid.dx);
A = obtain_linear_propagator(t_mask, f0, medium.sound_speed, kgrid.dx, set_current_A, 'mask_active', mask_active);

% Solve inverse problem
tic
if only_focus_opt
    % Preserve sonicated points only
    obs_ids = find(b_mask);
    ip.A = A(obs_ids, :);
    b_ip_des = b_des;

    ip.sq_beta = 0;%0.1; % constraint scaling factor
else
    % Take entire observation grid into account
    ip.A = A;
    b_ip_des = b_des_pl;
    ip.sq_beta = 0;%100; % constraint scaling factor
end

ip.p = pinv(ip.A) * b_ip_des; % Initial solution cosidering phases

ip.u = max(abs(ip.p)) * ones(size(ip.p)); % universal transducer amplitude

% Apply constraints (same transducer amplitude among elements)
ip.A = [ip.A; ip.sq_beta * eye(length(ip.p))];
b_ip_des = [b_ip_des; ip.sq_beta * ip.u];

% Solve phase retrieval problem
opts = struct;
opts.initMethod = 'custom';
opts.customx0 = ip.p;
[ip.p, outs, opts] = solvePhaseRetrieval(ip.A, ip.A', b_ip_des, [], opts);

ip.t_solve = toc;

% Evaluate obtained phase terms in forward simulation
ip.b = sim_exe(kgrid, medium, sensor, f0, ip.p, t_mask, sensor_mask, true, input_args, 'karray_t', karray_t);

%% Save Results in mat-file
if save_results
    current_datetime = string(datestr(now, 'yyyymmddHHMMSS'));
    res_filename = "results";
    if ~only_focus_opt
        ip.A = []; % A might be very large...
    end
    save(fullfile("Results", current_datetime + "_" + res_filename + ".mat"), ...
        "f0", "kgrid", "b_mask", "t_mask", "mask_active", "mask2el_delayFiles", "t_pos", "t_rot", "tr", "ip", "amp_in", "point_pos", "point_pos_m", ...
        "only_focus_opt", "input_args");
end

%% TR Results
if do_time_reversal
    plot_results(kgrid, tr.p, tr.b, t_pos, 'Time Reversal');
end

%% IP Results
if kgrid.dim == 2
    varargin = {};
else
    varargin = {'z_coord', point_pos.z(1)};
end

plot_results(kgrid, ip.p, ip.b, t_pos, 'Inverse Problem', varargin{:});

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
            b_tr_points = [b_tr_points, tr.b(point_pos.x(point), point_pos.y(point))];
            b_ip_points = [b_ip_points, ip.b(point_pos.x(point), point_pos.y(point))];
        end
    else
        for point = 1:length(point_pos.x)
            b_tr_points = [b_tr_points, tr.b(point_pos.x(point), point_pos.y(point), point_pos.z(point))];
            b_ip_points = [b_ip_points, ip.b(point_pos.x(point), point_pos.y(point), point_pos.z(point))];
        end
    end
    
    fprintf("\nInput Amplitudes (kPa):\n")
    disp(amp_in' * 1e-3)
    fprintf("\nTime Reversal Total Amplitudes (kPa):\n")
    disp(abs(b_tr_points) * 1e-3)
    fprintf("\nInverse Problem Total Amplitudes (kPa):\n")
    disp(abs(b_ip_points) * 1e-3)
    fprintf("\nTime Reversal Phase Angles (deg):\n")
    disp(angle(b_tr_points) / pi * 180)
    fprintf("\nInverse Problem Phase Angles (deg):\n")
    disp(angle(b_ip_points) / pi * 180)
end

