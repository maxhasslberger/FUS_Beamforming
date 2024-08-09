function plot_focalpt_waveform(f0,A_cells,p,init_ids)
    n = size(A_cells{1},1);
    nfreq = length(f0);
    T_m = find_mixed_period(f0);
    % Computational load can be managed by tweaking sampling freq and/or time vector length
    Fs = 10e6; 
    dt = 1/Fs;
    t = 0:dt:T_m;  
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

    % Plot waveform at the focal point
    focus_signals = signals(init_ids, :);
    figure;
    plot(t, focus_signals(1,:));
    title('Signal at Focal Point 1');
    xlabel('Time (s)');
    ylabel('Amplitude (Pa)');
    grid on;  

    figure;
    plot(t, focus_signals(2,:));
    title('Signal at Focal Point 2');
    xlabel('Time (s)');
    ylabel('Amplitude (Pa)');
    grid on;  

    % Plot the fft of the signal
    N = length(focus_signals(1,:));
    Y = fft(focus_signals(1,:));
    P2 = abs(Y/N);
    P1 = Y(1:(N/2)+1); % Plus 1 for Nyquist freq
    P1(2:end-1) = 2*P1(2:end-1); % Energy of negative spectrum is folded onto the positive spectrum
    k = 0:N/2; % One-sided fft
    freq = k * Fs/N;
    
    figure;
    stem(freq,P1)
    title("One-sided absolute value of focus signal fft");
    xlabel("Frequency (Hz)");
    ylabel("Amplitude (Pa)");
    xlim([0 600000])


end

function T_m = find_mixed_period(f0)
    T_array = 1 ./ f0;
    T_max = max(T_array);
    n = 1;
    found_soln = false;
    while ~found_soln
        if (sum(mod(n*T_max,T_array)) > 0.0)
            n = n+1;
        else
            T_m = n*T_max;
            found_soln = true;
        end
    end
end