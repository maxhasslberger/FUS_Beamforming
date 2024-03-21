clear;
close all;

%% Init

f0 = 500e3; % Hz - transducer frequency
n_dim = 2;
dx_factor = 1;
[kgrid, medium, ppp] = init_grid_medium(f0, 'n_dim', n_dim, 'dx_factor', 1 / dx_factor);
sensor = init_sensor(kgrid, ppp);
sensor_plane = ones(kgrid.Nx, kgrid.Ny);

only_focus_opt = true; % Optimize only focal spots or entire grid

%% Define Transducer Geometry

% karray_t = kWaveArray();

el1_offset = round(0.05 * kgrid.Nx); % grid points
el2_offset = round(0.05 * kgrid.Ny); % grid points
shift = 0; % grid points -> tangential shift

% Linear array
if only_focus_opt
    num_elements = 84;
    spacing = ceil(1 * dx_factor); % grid points between elements
    t_mask = create_linear_array(kgrid, num_elements, el1_offset, shift, spacing, false);
    % t_mask = t_mask + create_linear_array(kgrid, num_elements, el2_offset, shift, spacing, true); % Second (orthogonal) linear array
    % t_mask = t_mask + create_linear_array(kgrid, num_elements, round(kgrid.Nx - el1_offset), shift, spacing, false); % Second (antiparallel) linear array
    % t_mask = t_mask + create_linear_array(kgrid, num_elements, round(kgrid.Nx - el2_offset), shift, spacing, true); % Third (orthogonal) linear array
else
    num_elements = 42;
    spacing = ceil(2 * dx_factor); % grid points between elements
    t_mask = create_linear_array(kgrid, num_elements, el1_offset, shift, spacing, false);
    t_mask = t_mask + create_linear_array(kgrid, num_elements, el2_offset, shift, spacing, true); % Second (orthogonal) linear array
    % t_mask = t_mask + create_linear_array(kgrid, num_elements, round(kgrid.Nx - el1_offset), shift, spacing, false); % Second (antiparallel) linear array
    % t_mask = t_mask + create_linear_array(kgrid, num_elements, round(kgrid.Nx - el2_offset), shift, spacing, true); % Third (orthogonal) linear array
end

% % Curvature
% grid_size = [kgrid.Nx, kgrid.Ny];
% radius = kgrid.Nx; % grid points -> curvature
% diameter = round(0.6 * kgrid.Ny + 1 - mod(kgrid.Ny, 2)); % grid points -> one end to the other
% focus_pos = round(grid_size / 2);
% 
% t_mask = makeArc(grid_size, [el1_offset, round(kgrid.Ny / 2 + shift)], radius, diameter, focus_pos);
% t_mask = t_mask + makeArc(grid_size, [round(kgrid.Nx / 2 + shift), el2_offset], radius, diameter, focus_pos); % Second (orthogonal) curved array
% % t_mask = t_mask + makeArc(grid_size, [round(kgrid.Nx - el1_offset), round(kgrid.Ny / 2 + shift)], radius, diameter, focus_pos); % Second (antiparallel) curved array

% TODO: Define both transducer geometries in one (linear -> radius = inf;)

t_mask = t_mask > 0; % Return to logical in case of overlaps
% imagesc(t_mask, [-1 1])
% colormap(getColorMap);

%% Define (intracranial) Beamforming Pattern

if kgrid.dim == 2
    
    if ~only_focus_opt

        % Ring
        b_mask = makeDisc(kgrid.Nx, kgrid.Ny, kgrid.Nx/2, kgrid.Ny/2, round(0.2 * kgrid.Nx), false) ...
            - makeDisc(kgrid.Nx, kgrid.Ny, kgrid.Nx/2, kgrid.Ny/2, round(0.15 * kgrid.Nx), false);
        amp = 30000 * ones(sum(b_mask(:)), 1);
    else

        % Points
        point_posx = round([0.5, 0.2, 0.8] * kgrid.Nx);
        point_posy = round([0.6, 0.5, 0.2] * kgrid.Ny);
        amp = [10, 10, 10]' * 1e3;
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
else
    error("Not supported at the moment")
end


phase = zeros(length(amp), 1); % Zero phase for entire observation plane

b_des = amp .* exp(1j*phase); % only observed elements
b_des_pl = zeros(kgrid.Nx * kgrid.Ny, 1); % entire plane
b_des_pl(find(b_mask)) = b_des;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imagesc(b_mask + t_mask, [-1 1])
colormap(getColorMap);

% set simulation input options
input_args = {'PMLSize', 'auto', 'PMLInside', false, 'PlotPML', true, 'DisplayMask', b_mask + t_mask, 'RecordMovie', false};

%% Time Reversal

p_tr = sim_exe(kgrid, medium, sensor, f0, b_des, b_mask, t_mask, false, input_args);
p_tr = conj(p_tr);
b_tr = sim_exe(kgrid, medium, sensor, f0, p_tr, t_mask, sensor_plane, true, input_args);

%% Inverse Problem

% Obtain propagation operator
% A = linearPropagator_vs_acousticFieldPropagator(t_mask, f0, medium.sound_speed, kgrid.dx);
A = obtain_linear_propagator(t_mask, b_mask, f0, medium.sound_speed, kgrid.dx, only_focus_opt); % acousticFieldPropagator (Green's functions)

% Solve inverse problem
tic
if only_focus_opt
    p_ip = pinv(A) * b_des;
else

    p_ip = pinv(A) * b_des_pl;
    
    opts = struct;
    opts.initMethod = 'custom';
    opts.customx0 = p_ip;
    [p_ip, outs, opts] = solvePhaseRetrieval(A, A', b_des_pl, [], opts);
end

t_solve = toc;

b_ip = sim_exe(kgrid, medium, sensor, f0, p_ip, t_mask, sensor_plane, true, input_args);


%% Results
plot_results(kgrid, p_tr, b_tr, 'Time Reversal');
plot_results(kgrid, p_ip, b_ip, 'Inverse Problem');

disp("Time until solver converged: " + string(t_solve) + " s")

if only_focus_opt
    b_tr_points = [];
    b_ip_points = [];
    for point = 1:length(point_posx)
        b_tr_points = [b_tr_points, b_tr(point_posx(point), point_posy(point))];
        b_ip_points = [b_ip_points, b_ip(point_posx(point), point_posy(point))];
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

