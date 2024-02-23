
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
x_offset = 25;
spacing = 1;
t_mask = zeros(kgrid.Nx, kgrid.Ny);
start_index = kgrid.Ny/2 - round(num_elements/2) * spacing + 1;
t_mask(x_offset, start_index:spacing:start_index + num_elements * spacing - 1) = 1;


% if kgrid.dim == 2
%     n_elements = 10;
%     element_pitch = 10e-3; % m
%     translation = [0, 0];
%     rotation = 0;
%     element_width = 1e-3;
%     
%     % add rectangular array elements
%     for ind = 1:n_elements
%         
%         % set element y position
%         x_pos = 0 - (n_elements * element_pitch / 2 - element_pitch / 2) + (ind - 1) * element_pitch;
%         
%         % add element (set rotation angle to match the global rotation angle)
%         karray_t.addRectElement([x_pos, kgrid.y_vec(1)], element_width, element_width, rotation);
%         
%     end
% end
% 
% karray_t.setArrayPosition(translation, rotation)
% t_mask = karray_t.getArrayBinaryMask(kgrid);

%% Define (intracranial) Beamforming Pattern


if kgrid.dim == 2
    % Ring

    b_mask = makeDisc(kgrid.Nx, kgrid.Ny, kgrid.Nx/2, kgrid.Ny/2, 10, false) ...
        - makeDisc(kgrid.Nx, kgrid.Ny, kgrid.Nx/2, kgrid.Ny/2, 5, false);

%     karray_b = kWaveArray();
%     r = 25e-3;
%     dia = 50e-3;
%     pos = [0, r];
%     foc = [0, 1e-3];
% 
%     karray_b.addArcElement(+pos, r, dia, +foc)
%     karray_b.addArcElement(-pos, r, dia, -foc)
else
    
end

% b_mask = karray_b.getArrayBinaryMask(kgrid);
imagesc(b_mask)

amp = 30000 * ones(sum(b_mask(:)), 1);
phase = zeros(length(amp), 1); % Zero phase for entire observation plane

b_des = amp .* exp(1j*phase);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imagesc(b_mask + t_mask)

%% Time Reversal

p_tr = sim_exe(kgrid, medium, f0, karray_b, b_des, b_mask, t_mask, input_args);
p_tr = conj(p_tr);
b_tr = sim_exe(kgrid, medium, f0, karray_t, p_tr, t_mask, b_mask, input_args);

%% Inverse Problem

% Obtain propagation operator
A = obtain_propagation_operator(t_mask, b_mask); % acousticFieldPropagator (Green's functions) vs. angularSpectrum vs. focus function?

p_ip = pinv(A) * b_des;
b_ip = sim_exe(kgrid, medium, f0, p_ip, t_mask, b_mask, input_args);

%% Plot
plot_results(kgrid, b_tr, 'Time Reversal');
plot_results(kgrid, b_ip, 'Inverse Problem');


