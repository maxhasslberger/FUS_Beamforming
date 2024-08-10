function plot_focalpt_waveform(f0,A_cells,p,init_ids)
    n = size(A_cells{1},1);
    nfreq = length(f0);
    T_m = find_mixed_period(f0);
    % Computational load can be managed by tweaking sampling freq and/or time vector length
    Fs = 100e6; 
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
    % Take only focal point values
    amp_init = amp(init_ids, :, :);
    phi_init = phi(init_ids, :, :);
    % Broadcast w to match amp and phi
    w = reshape(w, 1, 1, nfreq);
    % Construct signals
    signals_tmp = amp_init .* cos(w .* t + phi_init);
    % Sum signals of different frequencies along the third dimension to get summed signals matrix
    signals = sum(signals_tmp, 3);

    % Plot waveform at the focal point
    figure;
    plot(t, signals(1,:));
    title('Signal at Focal Point 1');
    xlabel('Time (s)');
    ylabel('Amplitude (Pa)');
    grid on;  

    figure;
    plot(t, signals(2,:));
    title('Signal at Focal Point 2');
    xlabel('Time (s)');
    ylabel('Amplitude (Pa)');
    grid on;  

    % Plot the fft of the signal
    focus_signal_1 = signals(1,:);
    N = length(focus_signal_1);
    Y = fft(focus_signal_1);
    P2 = abs(Y)/N;
    P1 = P2(1:N/2+1); % Plus 1 for Nyquist freq
    P1(2:end-1) = 2*P1(2:end-1); % Energy of negative spectrum is folded onto the positive spectrum
    k = 0:N/2; % One-sided fft
    freq = k * Fs/N;
    phase = angle(Y(1:N/2+1));
    
    figure;
    subplot(2,1,1)
    stem(freq,P1)
    title("One-sided Amplitude Spectrum of Focus Signal 1");  
    ylabel("Amplitude (Pa)");
    xlim([300000 600000]);
    subplot(2,1,2);
    stem(freq,phase);
    title("One-sided Phase Spectrum of Focus Signal 1");
    xlabel('Frequency(Hz)');
    ylabel('Phase (rad)')
    xlim([300000 600000]);
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