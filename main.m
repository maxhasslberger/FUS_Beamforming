
%% Init
f0 = 500e3; % Hz - transducer frequency
[kgrid, medium, ppp] = init_grid_medium(f0, 'n_dim', 2);
sensor = init_sensor(kgrid, ppp);

% set simulation input options
input_args = {'PMLSize', 'auto', 'PMLInside', false, 'PlotPML', true, 'DisplayMask', 'off'};

%% Define Transducer Geometry
karray = kWaveArray('BLITolerance', 0.01, 'UpsamplingRate', 16);

% Annular
% karray.addAnnularArray(bowl_pos, source_roc, diameters, focus_pos);

% Linear Array
x_line = kgrid.x_vec(round(kgrid.Nx / 4));

if kgrid.dim == 2
    start_point = [x_line, kgrid.y_vec(round(kgrid.Ny / 4))];
    end_point = [x_line, kgrid.y_vec(round(kgrid.Ny * 3/4))];
else
    z_line = kgrid.z_vec(round(kgrid.Nz/ 2));
    start_point = [x_line, kgrid.y_vec(round(kgrid.Ny / 4)), z_line];
    end_point = [x_line, kgrid.y_vec(round(kgrid.Ny * 3/4)), z_line];
end
karray.addLineElement(start_point, end_point);

transducer_mask = karray.getArrayBinaryMask(kgrid);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Define (intracranial) Beamforming Pattern
if kgrid.dim == 2
    % Ring
    b_mask = makeDisc(kgrid.Nx, kgrid.Ny, kgrid.Nx/2, kgrid.Ny/2, 20, false) ...
        - makeDisc(kgrid.Nx, kgrid.Ny, kgrid.Nx/2, kgrid.Ny/2, 15, false);
else
    
end

%% Initialize Acoustic Signal
% amp = 
% phase =
p = amp * exp(1j*phase);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Source
source.p_mask = transducer_mask;
source_sig = createCWSignals(kgrid.t_array, f0, abs(p), angle(p));
source.p = karray.getDistributedSourceSignal(kgrid, source_sig);

%% Sensor
sensor.mask = b_mask;

%% Run Acoustic Simulation
if kgrid.dim == 2
    sensor_data = kspaceFirstOrder2DC(kgrid, medium, source, sensor, input_args{:});
else
    sensor_data = kspaceFirstOrder3DC(kgrid, medium, source, sensor, input_args{:});
end

p = post_processing(kgrid, f0, sensor_data.p);



plot_results(kgrid, p);
