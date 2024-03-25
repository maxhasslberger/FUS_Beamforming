clear;
close all;

%% Init

f0 = 500e3; % Hz - transducer frequency
n_dim = 2;
dx_factor = 1;
if n_dim == 2
    grid_size = [100, 100] * 1e-3; % m in [x, y] respectively
else
    grid_size = [150, 150, 100] * 1e-3; % m in [x, y, z] respectively
end
[kgrid, medium, ppp] = init_grid_medium(f0, grid_size, 'n_dim', n_dim, 'dx_factor', 1 / dx_factor);
[sensor, sensor_mask] = init_sensor(kgrid, ppp);

only_focus_opt = true; % Optimize only focal spots or entire grid
set_current_A = false; % Use precomputed propagation matrix - can be logical or a string containing the file name in Lin_Prop_Matrices

%% Define Transducer Geometry

if kgrid.dim == 2
    % Linear array
    t1_pos = -0.045; % m
    t2_pos = -0.045; % m

    el1_offset = round((t1_pos - kgrid.x_vec(1)) / kgrid.dx); % grid points
    el2_offset = round((t2_pos - kgrid.y_vec(1)) / kgrid.dy); % grid points
    if only_focus_opt
        num_elements = 50;
        shift = round(0e-3 / kgrid.dx); % m -> tangential shift in grid points
        spacing = ceil(1e-3 / kgrid.dx * dx_factor); % m -> grid points between elements
        t_mask = create_linear_array(kgrid, num_elements, el1_offset, shift, spacing, false);
        % t_mask = t_mask + create_linear_array(kgrid, num_elements, el2_offset, shift, spacing, true); % Second (orthogonal) linear array
        % t_mask = t_mask + create_linear_array(kgrid, num_elements, round(kgrid.Nx - el1_offset), shift, spacing, false); % Second (antiparallel) linear array
        % t_mask = t_mask + create_linear_array(kgrid, num_elements, round(kgrid.Nx - el2_offset), shift, spacing, true); % Third (orthogonal) linear array
    else
        num_elements = 25;
        shift1 = round(20e-3 / kgrid.dx); % m -> tangential shift in grid points
        shift2 = -shift1;
        spacing = ceil(2e-3 / kgrid.dx * dx_factor); % m -> grid points between elements
        t_mask = create_linear_array(kgrid, num_elements, el1_offset, shift1, spacing, false);
        t_mask = t_mask + create_linear_array(kgrid, num_elements, el2_offset, shift2, spacing, true); % Second (orthogonal) linear array
        % t_mask = t_mask + create_linear_array(kgrid, num_elements, round(kgrid.Nx - el1_offset), shift1, spacing, false); % Second (antiparallel) linear array
        % t_mask = t_mask + create_linear_array(kgrid, num_elements, round(kgrid.Nx - el2_offset), shift2, spacing, true); % Third (orthogonal) linear array
    end
    
    karray_t = [];
    t_mask = t_mask > 0; % Return to logical in case of overlaps

%     imagesc(t_mask, [-1 1])
%     colormap(getColorMap);
else
    % Planar Array
    t_name = "std";
    t1_pos = [-70, 20, 0]' * 1e-3; % m
    t1_rot = [0, 90, 0]'; % deg
    t2_pos = [20, -70, 0]' * 1e-3; % m
    t2_rot = [90, 0, 0]'; % deg

    t_pos = [t1_pos, t2_pos];
    t_rot = [t1_rot, t2_rot];

    karray_t = create_transducer(kgrid, t_name, t_pos, t_rot);
    t_mask = karray_t.getArrayBinaryMask(kgrid);

%     voxelPlot(double(t_mask))
end

%% Define (intracranial) Beamforming Pattern

if kgrid.dim == 2
    
    if ~only_focus_opt

        % Ring
        b_mask = makeDisc(kgrid.Nx, kgrid.Ny, round(0.7 * kgrid.Nx), round(0.7 * kgrid.Ny), round(0.2 * kgrid.Nx), false) ...
            - makeDisc(kgrid.Nx, kgrid.Ny, round(0.7 * kgrid.Nx), round(0.7 * kgrid.Ny), round(0.15 * kgrid.Nx), false);
        amp = 30000 * ones(sum(b_mask(:)), 1);
    else

        % Points
        point_posx = round([0.2, 0.5, 0.6] * kgrid.Nx);
        point_posy = round([0.5, 0.6, 0.3] * kgrid.Ny);
        point_posz = [];
        amp = [0, 100, 200]' * 1e3;
%         amp = 30000 * ones(length(point_posx), 1);
    
        % Assign amplitude acc. to closest position
        idx = sub2ind([kgrid.Nx, kgrid.Ny], point_posx, point_posy);
        [~, order] = sort(idx);
        amp = amp(order);
    
        b_mask = zeros(kgrid.Nx, kgrid.Ny);

        for point = 1:length(point_posx)
            b_mask(point_posx(point), point_posy(point)) = 1;
        end
    end

    imagesc(b_mask + t_mask, [-1 1])
    colormap(getColorMap);
else
    only_focus_opt = true;

    % Points
    point_posx = round([0.2, 0.7, 0.8] * kgrid.Nx);
    point_posy = round([0.5, 0.8, 0.5] * kgrid.Ny);
    point_posz = round([0.5, 0.5, 0.5] * kgrid.Nz);
    amp = [10, 10, 10]' * 1e3;

    % Assign amplitude acc. to closest position
    idx = sub2ind([kgrid.Nx, kgrid.Ny, kgrid.Nz], point_posx, point_posy, point_posz);
    [~, order] = sort(idx);
    amp = amp(order);

    b_mask = zeros(kgrid.Nx, kgrid.Ny, kgrid.Nz);

    for point = 1:length(point_posx)
        b_mask(point_posx(point), point_posy(point), point_posz(point)) = 1;
    end

    voxelPlot(double(t_mask + b_mask))
end

% Create desired signal
phase = zeros(length(amp), 1); % Zero phase for entire observation plane

b_des = amp .* exp(1j*phase); % only observed elements

b_des_pl = zeros(kgrid.Nx * kgrid.Ny, 1); % entire plane
b_des_pl(find(b_mask)) = b_des;

% set simulation input options
input_args = {'PMLSize', 'auto', 'PMLInside', false, 'PlotPML', true, 'DisplayMask', b_mask + t_mask, 'RecordMovie', false};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Time Reversal

p_tr = sim_exe(kgrid, medium, sensor, f0, b_des, b_mask, t_mask, false, input_args);
p_tr = max(abs(p_tr)) * exp(-1j * angle(p_tr)); % All elements with same amplitude
b_tr = sim_exe(kgrid, medium, sensor, f0, p_tr, t_mask, sensor_mask, true, input_args, 'karray_t', karray_t);

%% Inverse Problem

% Obtain propagation operator
% A = linearPropagator_vs_acousticFieldPropagator(t_mask, f0, medium.sound_speed, kgrid.dx);
A = obtain_linear_propagator(t_mask, b_mask, f0, medium.sound_speed, kgrid.dx, only_focus_opt, set_current_A); % -> acousticFieldPropagator (Green's functions)

% Solve inverse problem
tic
if only_focus_opt
    b_ip_des = b_des;
    sq_beta = 0;%0.1; % constraint scaling factor
else
    b_ip_des = b_des_pl;
    sq_beta = 0;%100; % constraint scaling factor
end

p_ip = pinv(A) * b_ip_des; % Initial solution cosidering phases

u = max(abs(p_ip)) * ones(size(p_ip)); % universal transducer amplitude

% Apply constraints (same transducer amplitude among elements)
A = [A; sq_beta * eye(length(p_ip))];
b_ip_des = [b_ip_des; sq_beta * u];

% Solve phase retrieval problem
opts = struct;
opts.initMethod = 'custom';
opts.customx0 = p_ip;
[p_ip, outs, opts] = solvePhaseRetrieval(A, A', b_ip_des, [], opts);

t_solve = toc;

% Evaluate obtained phase terms in forward simulation
b_ip = sim_exe(kgrid, medium, sensor, f0, p_ip, t_mask, sensor_mask, true, input_args, 'karray_t', karray_t);

%% Save Results in mat-file
% f0
% kgrid
% b_mask
% t_mask
% t1_pos
% t2_pos

tr.p_tr = p_tr;
tr.b_tr = b_tr;

ip.A = A;
ip.p_ip = p_ip;
ip.b_ip = b_ip;
ip.t_solve = t_solve;
ip.u = u;
ip.sq_beta = sq_beta;

point_pos.x = point_posx;
point_pos.y = point_posy;
point_pos.z = point_posz;

%% Results
plot_results(kgrid, p_tr, b_tr, 'Time Reversal');
plot_results(kgrid, p_ip, b_ip, 'Inverse Problem');

% Metrics evaluation
disp("Time until solver converged: " + string(t_solve) + " s")
fprintf("\nDesired Transducer Amplitude (kPa):\n")
disp(u(1) * 1e-3)
fprintf("\nTransducer Amplitude mean deviation (Pa):\n")
disp(mean(abs(abs(p_ip) - u)))

% Evaluation of each defined point
if only_focus_opt
    b_tr_points = [];
    b_ip_points = [];

    if kgrid.dim == 2
        for point = 1:length(point_posx)
            b_tr_points = [b_tr_points, b_tr(point_posx(point), point_posy(point))];
            b_ip_points = [b_ip_points, b_ip(point_posx(point), point_posy(point))];
        end
    else
        for point = 1:length(point_posx)
            b_tr_points = [b_tr_points, b_tr(point_posx(point), point_posy(point), point_posz(point))];
            b_ip_points = [b_ip_points, b_ip(point_posx(point), point_posy(point), point_posz(point))];
        end
    end
    
    fprintf("\nInput Amplitudes (kPa):\n")
    disp(amp' * 1e-3)
    fprintf("\nTime Reversal Total Amplitudes (kPa):\n")
    disp(abs(b_tr_points) * 1e-3)
    fprintf("\nInverse Problem Total Amplitudes (kPa):\n")
    disp(abs(b_ip_points) * 1e-3)
    fprintf("\nTime Reversal Phase Angles (kPa):\n")
    disp(angle(b_tr_points) * 1e-3)
    fprintf("\nInverse Problem Phase Angles (kPa):\n")
    disp(angle(b_ip_points) * 1e-3)
end

