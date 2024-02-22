function sensor = init_sensor(kgrid, ppp)

record_periods = 3;

% record the pressure
sensor.record = {'p'};

% average only the final few periods when the field is in steady state
sensor.record_start_index = kgrid.Nt - record_periods * ppp + 1;

end
