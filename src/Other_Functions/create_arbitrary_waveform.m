% clear;
% close all;
 
use_cw_signals = false;


% -------------------------------------------------------
% >>>>> Option 1: Constructing signals from scratch <<<<<
% -------------------------------------------------------

if ~use_cw_signals
    f0 = [450:10:490] * 1e3; % ctr freq = 470 kHz
    nfreq = length(f0);
    T = 1/(min(f0));
    % Computational load can be managed by tweaking sampling freq and/or time vector length
    Fs = 100e6; 
    dt = 1/Fs;
    t = 0:dt:100*T;  
    w = 2 * pi * f0; % Angular frequencies
    amp = 100 * 1e3; % Pa
    
    % Construct signals
    signals_tmp = amp .* cos(w' .* t);
    % Sum signals of different frequencies along the third dimension to get summed signals matrix
    result_signal = sum(signals_tmp, 1);
    
    figure;
    plot(t, result_signal);
    title('Arbitrary Signal');
    xlabel('Time (s)');
    ylabel('Amplitude');
    grid on;

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