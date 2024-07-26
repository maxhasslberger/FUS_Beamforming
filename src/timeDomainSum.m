function y = timeDomainSum(f0,A_cells,p)
    n = size(A_cells{1},1);
    nfreq = length(f0);
    T = 1/(min(f0));
    % Computational load can be managed by tweaking sampling freq and/or time vector length
    Fs = 10e6; 
    dt = 1/Fs;
    t = 0:dt:2*T;  
    w = 2 * pi * f0; % Angular frequencies

    % Find amps and phases
    amp = zeros(n,1,nfreq);
    phi = zeros(n,1,nfreq);
    for i = 1:nfreq
        y_temp = A_cells{i} * p(:, i);
        amp(:,1,i) = abs(y_temp); 
        phi(:,1,i) = angle(y_temp);
    end

    % Broadcast w to match amp and phi
    w = reshape(w, 1, 1, nfreq);
    % Construct signals
    signals_tmp = amp .* cos(w .* t + phi);
    % Sum signals of different frequencies along the third dimension to get summed signals matrix
    signals = sum(signals_tmp, 3);

    % Compute amplitudes
    y = (max(signals, [], 2) - min(signals, [], 2)) / 2;
end
