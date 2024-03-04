function p = sim_exe(kgrid, medium, f0, p_in, source_mask, sensor_mask, final_sim, input_args)

%% Source

source.p_mask = source_mask;
source.p = createCWSignals(kgrid.t_array, f0, abs(p_in), angle(p_in));
% source.p = karray.getDistributedSourceSignal(kgrid, source_sig);

%% Sensor

sensor.mask = sensor_mask;

%% Run Acoustic Simulation

if kgrid.dim == 2
    sensor_data = kspaceFirstOrder2DC(kgrid, medium, source, sensor, input_args{:});
else
    sensor_data = kspaceFirstOrder3DC(kgrid, medium, source, sensor, input_args{:});
end

p = post_processing(kgrid, f0, sensor_data, final_sim);

end
