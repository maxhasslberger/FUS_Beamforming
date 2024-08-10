function [amp_in, phase_in] = get_amp_phase_mask(kgrid, f0, p, t_mask, karray_t)

if ~isempty(karray_t)
    % Obtain karray signals in TD
    source_signal = createCWSignals(kgrid.t_array, f0, abs(p), angle(p));
    t_mask = karray_t.getArrayBinaryMask(kgrid);
    signal_td = karray_t.getDistributedSourceSignal(kgrid, source_signal);

    % Convert to complex excitations
    [amp, phase] = extractAmpPhase(signal_td, 1/kgrid.dt, f0, ...
    'Dim', 2, 'Window', 'Rectangular', 'FFTPadding', 1);
    p = amp .* exp(1j*phase);
end

% Assign excitation in space
signal = double(t_mask);
signal(t_mask) = p;
amp_in = abs(signal);
phase_in = angle(signal);

end

