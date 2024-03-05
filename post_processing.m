function p = post_processing(kgrid, f0, data, final_sim)

% extract amplitude and phase from the sensor data
[amp, phase] = extractAmpPhase(data, 1/kgrid.dt, f0, ...
    'Dim', 2, 'Window', 'Rectangular', 'FFTPadding', 1);

% reshape data
Nx = kgrid.Nx;
Ny = kgrid.Ny;

if final_sim

    if kgrid.dim == 2
        amp = reshape(amp, Nx, Ny);
        phase = reshape(phase, Nx, Ny);
    else % kgrid.dim == 3
        Nz = kgrid.Nz;
    
        amp = reshape(amp, Nx, Ny, Nz);
        phase = reshape(phase, Nx, Ny, Nz);
    end
else
%     amp = max(amp) * ones(length(amp), 1); % Same amplitude for each element
end

p = amp .* exp(1j*phase);

end
