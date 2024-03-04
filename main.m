clear;
close all;

%% Init

f0 = 500e3; % Hz - transducer frequency
n_dim = 2;
dx_factor = 1.0;
[kgrid, medium, ppp] = init_grid_medium(f0, 'n_dim', n_dim, 'dx_factor', 1 / dx_factor);
sensor = init_sensor(kgrid, ppp);

% set simulation input options
input_args = {'PMLSize', 'auto', 'PMLInside', false, 'PlotPML', true, 'DisplayMask', 'off'};

%% Define Transducer Geometry

% karray_t = kWaveArray();

% Annular
% karray.addAnnularArray(bowl_pos, source_roc, diameters, focus_pos);

% Linear Array
num_elements = 84 * dx_factor;
x_offset = 0.05 * kgrid.Ny; % grid points
spacing = 1; % grid points between elements

t_mask = create_linear_array(kgrid, num_elements, x_offset, spacing); % TODO: More elements, circular
sensor_plane = ones(kgrid.Nx, kgrid.Ny);

%% Define (intracranial) Beamforming Pattern


if kgrid.dim == 2

    % Ring
    b_mask = makeDisc(kgrid.Nx, kgrid.Ny, kgrid.Nx/2, kgrid.Ny/2, 20, false) ...
        - makeDisc(kgrid.Nx, kgrid.Ny, kgrid.Nx/2, kgrid.Ny/2, 15, false);

    % Point
    point_pos = [0.5 * kgrid.Nx, 0.6 * kgrid.Ny];
    point_pos2 = [0.2 * kgrid.Nx, 0.5 * kgrid.Ny];

    b_mask = zeros(kgrid.Nx, kgrid.Ny);
    b_mask(point_pos(1), point_pos(2)) = 1;
    b_mask(point_pos2(1), point_pos2(2)) = 1;
%     b_mask(point_pos(1), kgrid.Ny - point_pos(2)) = 1;
else
    error("Not supported at the moment")
end


amp = 30000 * ones(sum(b_mask(:)), 1);
phase = zeros(length(amp), 1); % Zero phase for entire observation plane

b_des = amp .* exp(1j*phase); % only observed elements
b_des_pl = zeros(kgrid.Nx * kgrid.Ny, 1); % entire plane
b_des_pl(find(b_mask)) = b_des;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imagesc(b_mask + t_mask)

%% Time Reversal

p_tr = sim_exe(kgrid, medium, f0, b_des, b_mask, t_mask, false, input_args);
p_tr = conj(p_tr);
b_tr = sim_exe(kgrid, medium, f0, p_tr, t_mask, sensor_plane, true, input_args);

%% Inverse Problem

% Obtain propagation operator
% A = linearPropagator_vs_acousticFieldPropagator(t_mask, f0, medium.sound_speed, kgrid.dx);
A = obtain_linear_propagator(t_mask, b_mask, f0, medium.sound_speed, kgrid.dx); % acousticFieldPropagator (Green's functions)

p_ip = pinv(A) * b_des;
b_ip = sim_exe(kgrid, medium, f0, p_ip, t_mask, sensor_plane, true, input_args);

%% Plot
plot_results(kgrid, p_tr, b_tr, 'Time Reversal');
plot_results(kgrid, p_ip, b_ip, 'Inverse Problem');

% b_tr(point_pos(1), point_pos(2))
% b_ip(point_pos(1), point_pos(2))


