
%% Init

f0 = 500e3; % Hz - transducer frequency
n_dim = 2;
[kgrid, medium, ppp] = init_grid_medium(f0, 'n_dim', n_dim);
sensor = init_sensor(kgrid, ppp);

% set simulation input options
input_args = {'PMLSize', 'auto', 'PMLInside', false, 'PlotPML', true, 'DisplayMask', 'off'};

%% Define Transducer Geometry

% karray_t = kWaveArray();

% Annular
% karray.addAnnularArray(bowl_pos, source_roc, diameters, focus_pos);

% Linear Array
num_elements = 21;
x_offset = 25; % grid points
spacing = 1; % grid points between elements

t_mask = create_linear_array(kgrid, num_elements, x_offset, spacing);


%% Define (intracranial) Beamforming Pattern


if kgrid.dim == 2

    % Ring
    b_mask = makeDisc(kgrid.Nx, kgrid.Ny, kgrid.Nx/2, kgrid.Ny/2, 10, false) ...
        - makeDisc(kgrid.Nx, kgrid.Ny, kgrid.Nx/2, kgrid.Ny/2, 5, false);
else
    
end


amp = 30000 * ones(sum(b_mask(:)), 1);
phase = zeros(length(amp), 1); % Zero phase for entire observation plane

b_des = amp .* exp(1j*phase);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imagesc(b_mask + t_mask)

%% Time Reversal

% p_tr = sim_exe(kgrid, medium, f0, b_des, b_mask, t_mask, input_args);
% p_tr = conj(p_tr);
% b_tr = sim_exe(kgrid, medium, f0, p_tr, t_mask, b_mask, input_args);

%% Inverse Problem

% Obtain propagation operator
A = obtain_linear_propagator(t_mask, b_mask, f0, medium.sound_speed, kgrid.dx); % acousticFieldPropagator (Green's functions)

p_ip = pinv(A) * b_des;
b_ip = sim_exe(kgrid, medium, f0, p_ip, t_mask, b_mask, input_args);

%% Plot
plot_results(kgrid, b_tr, 'Time Reversal');
plot_results(kgrid, b_ip, 'Inverse Problem');


