% clear;
% close all;
 
use_cw_signals = false; % use kWave continuous wave function
perform_fft = false; % perform fft on summed signal
smoothen_pulse = true; % smooth pulses

% -------------------------------------------------------
% >>>>> Option 1: Constructing signals from scratch <<<<<
% -------------------------------------------------------

if ~use_cw_signals
    f0 = [450:10:490] * 1e3; % ctr freq = 470 kHz
    nfreq = length(f0);
    T = 1/(min(f0));
    % Computational load can be managed by tweaking sampling freq and/or time vector length
    Fs = 10e6; 
    dt = 1/Fs;
    % t = 0:dt:60*T;  
    w = 2 * pi * f0; % Angular frequencies
    amp = 20 * 1e3; % Pa

    npulses = 3;
    pulse_train_duration = 50 * T;
    t = 0:dt:(pulse_train_duration);
    pulse_duration = pulse_train_duration / npulses;
    duty_cycle = 0.5;
    
    % Construct signals
    signals_tmp = amp .* cos(w' .* t);
    % Sum signals of different frequencies along the first dimension to get summed signals matrix
    result_signal = sum(signals_tmp, 1);
    
    figure;
    plot(t, result_signal);
    title('Arbitrary Signal');
    xlabel('Time (s)');
    ylabel('Amplitude (Pa)');
    grid on;

    envelope = zeros(1, length(t));
    for pulse = 0:(npulses-1)
        sonicated_duration = duty_cycle * pulse_duration;
        pulse_start = pulse * pulse_duration;
        sonic_pulse_mid = pulse_start + sonicated_duration/2;
        sonic_pulse_end = pulse_start + sonicated_duration;

        if smoothen_pulse
            % Rising edge
            t_rise = t((t >= pulse_start) & (t < sonic_pulse_mid)) - pulse_start;
            envelope((t >= pulse_start) & (t < sonic_pulse_mid)) = 0.5 * (1 - cos(pi * t_rise / (sonicated_duration/2)));
            
            % Falling edge
            t_fall = t((t >= sonic_pulse_mid) & (t < sonic_pulse_end)) - sonic_pulse_mid;
            envelope((t >= sonic_pulse_mid) & (t < sonic_pulse_end)) = 0.5 * (1 + cos(pi * t_fall / (sonicated_duration/2)));
        else
            envelope((t>=pulse_start) & (t < sonic_pulse_end)) = 1;
        end
    end
    smooth_result_signal = result_signal .* envelope;

    % plot envelope
    figure;
    plot(t,envelope);
    title('envelope')
    grid on;

    % Plot
    figure;
    plot(t, smooth_result_signal);
    title('Burst Signal');
    xlabel('Time (s)');
    ylabel('Amplitude');
    grid on;

    if perform_fft
        % Compute the FFT of the resulting signal
        n = length(result_signal);
        Y = fft(result_signal);
        P2 = abs(Y/n); % Two-sided spectrum
        P1 = P2(1:n/2+1); % Single-sided spectrum
        P1(2:end-1) = 2*P1(2:end-1); % Energy of negative spectrum is folded onto the positive spectrum
        f_fft = Fs*(0:(n/2))/n;
    
        % Find peaks in the FFT result
        [peaks, locs] = findpeaks(P1, 'MinPeakHeight', 0.01*max(P1));
    
        % Plot the amplitude spectrum
        figure;
        plot(f_fft, P1);
        hold on;
        plot(f_fft(locs), peaks, 'ro');
        title('Amplitude Spectrum');
        xlabel('Frequency (Hz)');
        ylabel('Amplitude');
        grid on;
    
        % Zoom in on the frequency range of interest
        xlim([400e3 550e3]);
    end

% ----------------------------------------------------------------
% >>>>> Option 2: Constructing signals using createCWSignals <<<<<
% ----------------------------------------------------------------

else

    % define sampling parameters
    f0 = [450:10:490] * 1e3; % ctr freq = 470 kHz
    T = 1/(min(f0));
    Fs = 100e6;
    dt = 1/Fs;
    t_array = 0:dt:100*T;
    
    % define amplitude and phase
    amp = 1;
    phase = 0;
    
    % create signals
    cw_signal = createCWSignals(t_array, f0(1), amp, phase);
    for i = 2:length(f0)
        cw_signal = cw_signal + createCWSignals(t_array, f0(i), amp, phase);
    end
    
    % Plot the resulting summed signal
    figure;
    plot(t_array, cw_signal);
    title('Summed Signal');
    xlabel('Time (s)');
    ylabel('Amplitude');
    grid on;
    
    % stackedPlot(cw_signal);
end

