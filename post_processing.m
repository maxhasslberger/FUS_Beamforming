function p = post_processing(kgrid, f0, data, final_sim)

% reshape data
Nx = kgrid.Nx;
Ny = kgrid.Ny;

if final_sim
    % Only max amplitude matters
    amp = max(abs(data), [], 2);
    phase = zeros(length(amp), 1);

    if kgrid.dim == 2
        amp = reshape(amp, Nx, Ny);
        phase = reshape(phase, Nx, Ny);
    else % kgrid.dim == 3
        Nz = kgrid.Nz;
    
        amp = reshape(amp, Nx, Ny, Nz);
        phase = reshape(phase, Nx, Ny, Nz);
    end
else
    % extract amplitude and phase from the sensor data
    [amp, phase] = extractAmpPhase(data, 1/kgrid.dt, f0, ...
        'Dim', 2, 'Window', 'Rectangular', 'FFTPadding', 1);
end

p = amp .* exp(1j*phase);

end
