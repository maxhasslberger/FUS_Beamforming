function [sensor, sensor_mask] = init_sensor(kgrid, ppp)

record_periods = 3;

% record the pressure
% sensor.record = {'p'};

% average only the final few periods when the field is in steady state
sensor.record_start_index = kgrid.Nt - record_periods * ppp + 1;

if kgrid.dim == 2
    sensor_mask = ones(kgrid.Nx, kgrid.Ny);
else
    sensor_mask = ones(kgrid.Nx, kgrid.Ny, kgrid.Nz);
end

end
