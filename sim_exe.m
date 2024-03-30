function p = sim_exe(kgrid, medium, sensor, f0, p_in, source_mask, sensor_mask, final_sim, input_args, varargin)

karray_t = [];

if ~isempty(varargin)
    for arg_idx = 1:2:length(varargin)
        switch varargin{arg_idx}
            case 'karray_t'
                karray_t = varargin{arg_idx+1};
            otherwise
                error('Unknown optional input.');
        end
    end
end

%% Source

source_signal = createCWSignals(kgrid.t_array, f0, abs(p_in), angle(p_in));
if isempty(karray_t)
    source.p_mask = source_mask;
    source.p = source_signal;
else
    source.p_mask = karray_t.getArrayBinaryMask(kgrid);
    source.p = karray_t.getDistributedSourceSignal(kgrid, source_signal);
%     source.p = source_signal;
end

%% Sensor

sensor.mask = sensor_mask;

%% Run Acoustic Simulation
% Linux (in binaries folder): chmod +x ./kspaceFirstOrder3D-OMP or chmod +x ./kspaceFirstOrder-OMP
if kgrid.dim == 2
    sensor_data = kspaceFirstOrder2DC(kgrid, medium, source, sensor, input_args{:});
else
    sensor_data = kspaceFirstOrder3DC(kgrid, medium, source, sensor, input_args{:});
end

p = post_processing(kgrid, f0, sensor_data, final_sim);

end
