function p = post_processing(kgrid, f0, data)

% extract amplitude and phase from the sensor data
[amp, phase] = extractAmpPhase(data, 1/kgrid.dt, f0, ...
    'Dim', 2, 'Window', 'Rectangular', 'FFTPadding', 1);

% reshape data
if kgrid.dim == 2
    amp = reshape(amp, Nx, Ny);
    phase = reshape(phase, Nx, Ny);
else % kgrid.dim == 3
    amp = reshape(amp, Nx, Ny, Nz);
    phase = reshape(phase, Nx, Ny, Nz);
end

p = amp * exp(1j*phase);

end
