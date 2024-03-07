clear;
close all;

%% Init

f0 = 500e3; % Hz - transducer frequency
n_dim = 2;
dx_factor = 1;
[kgrid, medium, ppp] = init_grid_medium(f0, 'n_dim', n_dim, 'dx_factor', 1 / dx_factor);
sensor = init_sensor(kgrid, ppp);
sensor_plane = ones(kgrid.Nx, kgrid.Ny);

% set simulation input options
input_args = {'PMLSize', 'auto', 'PMLInside', false, 'PlotPML', true, 'DisplayMask', 'off'};

%% Define Transducer Geometry

% karray_t = kWaveArray();

% Annular
% karray.addAnnularArray(bowl_pos, source_roc, diameters, focus_pos);

% Linear Array
num_elements = 84;
el1_offset = round(0.05 * kgrid.Nx); % grid points
el2_offset = round(0.05 * kgrid.Ny); % grid points
spacing = ceil(1 * dx_factor); % grid points between elements
shift = 0; % grid points

% Create Transducer mask
t_mask = create_linear_array(kgrid, num_elements, el1_offset, shift, spacing, false);
% t_mask = t_mask + create_linear_array(kgrid, num_elements, el2_offset, shift, spacing, true); % Second (orthogonal) linear array
% t_mask = t_mask + create_linear_array(kgrid, num_elements, round(kgrid.Nx - el1_offset), shift, spacing, false); % Second (antiparallel) linear array
t_mask = t_mask > 0; % Return to logical in case of overlaps

%% Define (intracranial) Beamforming Pattern

if kgrid.dim == 2

%     % Ring
%     b_mask = makeDisc(kgrid.Nx, kgrid.Ny, kgrid.Nx/2, kgrid.Ny/2, 20, false) ...
%         - makeDisc(kgrid.Nx, kgrid.Ny, kgrid.Nx/2, kgrid.Ny/2, 15, false);
%     amp = 30000 * ones(sum(b_mask(:)), 1);

    % Points
    point_posx = round([0.5, 0.2, 0.8] * kgrid.Nx);
    point_posy = round([0.6, 0.5, 0.2] * kgrid.Ny);
    amp = [10, 10, 10]' * 1e3;
%     amp = 30000 * ones(length(point_posx), 1);

    % Assign amplitude acc. to closest position
    idx = sub2ind([kgrid.Nx, kgrid.Ny], point_posx, point_posy);
    [~, order] = sort(idx);
    amp = amp(order);

    b_mask = zeros(kgrid.Nx, kgrid.Ny);
    for point = 1:length(point_posx)
        b_mask(point_posx(point), point_posy(point)) = 1;
    end
else
    error("Not supported at the moment")
end


phase = zeros(length(amp), 1); % Zero phase for entire observation plane

b_des = amp .* exp(1j*phase); % only observed elements
b_des_pl = zeros(kgrid.Nx * kgrid.Ny, 1); % entire plane
b_des_pl(find(b_mask)) = b_des;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imagesc(b_mask + t_mask)

%% Time Reversal

p_tr = sim_exe(kgrid, medium, sensor, f0, b_des, b_mask, t_mask, false, input_args);
p_tr = conj(p_tr);
b_tr = sim_exe(kgrid, medium, sensor, f0, p_tr, t_mask, sensor_plane, true, input_args);

%% Inverse Problem

% Obtain propagation operator
% A = linearPropagator_vs_acousticFieldPropagator(t_mask, f0, medium.sound_speed, kgrid.dx);
A = obtain_linear_propagator(t_mask, b_mask, f0, medium.sound_speed, kgrid.dx); % acousticFieldPropagator (Green's functions)

p_ip = pinv(A) * b_des;
b_ip = sim_exe(kgrid, medium, sensor, f0, p_ip, t_mask, sensor_plane, true, input_args);

%% Results
plot_results(kgrid, p_tr, b_tr, 'Time Reversal');
plot_results(kgrid, p_ip, b_ip, 'Inverse Problem');


b_tr_points = [];
b_ip_points = [];
for point = 1:length(point_posx)
    b_tr_points = [b_tr_points, b_tr(point_posx(point), point_posy(point))];
    b_ip_points = [b_ip_points, b_ip(point_posx(point), point_posy(point))];
end

amp
abs(b_tr_points)
abs(b_ip_points)
angle(b_tr_points)
angle(b_ip_points)


