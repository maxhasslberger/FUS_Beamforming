clear;
close all;

% 50 microsecond burst
% 150 microsecond padding


% Modes
compute_waveform = true;
plot_pulse_durations = false;
use_cw_signals = false;
use_multi_freq = false;

% Waveform options
perform_fft = false; % perform fft on summed signal
smoothen_pulse = false; % smooth pulses
zero_padding = false; % apply zero padding to beginning of the signal
apply_factoring = false; 
apply_envelope = false;
save_csv = true; % export signal to csv file

% ----------------------------
% Construct arbitrary waveform
% ----------------------------

if compute_waveform  
        % Frequencies
        ctr_freq = 500; % kHz
        if use_multi_freq
            lowest_freq = ctr_freq - 10;
            highest_freq = ctr_freq + 10;
            f0 = [lowest_freq:5:highest_freq] * 1e3; 
        else 
            f0 = ctr_freq * 1e3;
        end
        nfreq = length(f0);
        w = 2 * pi * f0; % Angular frequencies        
        % Period
        if use_multi_freq
            T_m = find_mixed_period(f0);
        end
        T = 1/(min(f0));
        % Computational load can be managed by tweaking sampling freq and/or time vector length
        Fs = 100e6;
        dt = 1/Fs;
        % t = 0:dt:60*T;          
        % Amplitudes and phases of each freq
        if use_multi_freq
            amp_2d = [20, 20, 20, 20, 20] * 1e3; % kPa
            phi_2d = [0,0.2,0.4,0.6,0.8];
        else
            amp = 20*1e3;
            phi = 0;
        end
        if apply_factoring
            factors = [1/.97, 1/.99, 1, 1/.99, 1/.97];
            amp_2d_factored = amp_2d .* factors; % Pa                         
            amp = reshape(amp_2d_factored, 1, 1, nfreq);                    
        else                  
            if use_multi_freq
                amp = reshape(amp_2d, 1, 1, nfreq);  
            end
        end
        % Reshape vars to stack each frequency along third dim
        if use_multi_freq
            w = reshape(w, 1, 1, nfreq);
            phi = reshape(phi_2d, 1, 1, nfreq);                
            pulse_train_duration = T_m;
        else
            % pulse_train_duration = 30*T;
            pulse_train_duration = 200* 1e-6;
        if zero_padding; padding_duration = 10 * T; else; padding_duration = 0; end;
        t = 0:dt:(pulse_train_duration + padding_duration);
        npulses = 1;
        pulse_duration = pulse_train_duration / npulses;
        duty_cycle = 1;
        
        % Construct signals
        if use_multi_freq
            signals_tmp = amp .* cos(w .* t + phi);        
            % Sum signals of different frequencies along the first dimension to get summed signals matrix
            result_signal = sum(signals_tmp, 3);
        else
            result_signal = amp * sin(w .* t + phi);
        end
        % if zero_padding
        %     samples_zero = length(0:dt:padding_duration);
        %     result_signal(1:samples_zero) = 0; 
        % end
        time_zero = 150 * 1e-6;
        % time_zero = 10*T;
        t_zero = pulse_train_duration +dt:dt:pulse_train_duration + time_zero;
        zero_signal = zeros(1,length(t_zero));
                
        result_signal = [result_signal, zero_signal];
        total_time = pulse_train_duration + time_zero;
        t = 0:dt:total_time;
        figure;
        plot(t, result_signal);
        title('Arbitrary Signal');
        xlabel('Time (s)');
        ylabel('Amplitude (Pa)');
        grid on;
        ylim([-1.5e5 1.5e5]);
        
    
        if apply_envelope
            envelope = zeros(1, length(t));
            for pulse = 0:(npulses-1)
                sonicated_duration = duty_cycle * pulse_duration;
                pulse_start = padding_duration + pulse * pulse_duration;
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
        
            if zero_padding
                % smooth_result_signal = [zeros(1,50) smooth_result_signal];        
            end
        
            % Plot envelope
            figure;
            plot(t,envelope);
            title('envelope')
            grid on;      
            
            % Plot burst
            figure;
            plot(t, smooth_result_signal);
            title('Burst Signal');
            xlabel('Time (s)');
            ylabel('Amplitude');
            grid on;           
        end      
    
        if save_csv
            time = t';
            value = result_signal';
            signal_data = table(time, value);
            
            writetable(signal_data, fullfile("../..","Transducer_Delay_Files","arbitrary_signal_singlefreq_500khz.csv"));
            if ~use_multi_freq
                amplitude = amp;
                params = table(f0, amplitude, phi);
                writetable(params,fullfile("../..","Transducer_Delay_Files","single_freq" + num2str(f0) + "_params.csv"));
            else
                amplitude = amp_2d';
                phase_shift = phi_2d';
                freq = f0';
                unfactored_params = table(freq,amplitude,phase_shift);
                writetable(unfactored_params, fullfile("../..","Transducer_Delay_Files","unfactored_signal_params.csv"));            

                if apply_factoring
                    factored_amplitude = amp_2d_factored';
                    factored_params = table(freq,factored_amplitude,phase_shift);
                    writetable(factored_params, fullfile("../..","Transducer_Delay_Files","factored_signal_params.csv"));
                end
            end
        end
    
    

        if perform_fft
            % Compute the FFT of the resulting signal
            n = length(result_signal);
            Y = fft(result_signal);
            P2 = abs(Y/n); % Two-sided spectrum
            P1 = P2(1:n/2+1); % Single-sided spectrum
            P1(2:end-1) = 2*P1(2:end-1); % Energy of negative spectrum is folded onto the positive spectrum
            f_fft = Fs*(0:(n/2))/n;
            phase = angle(Y(1:n/2+1));
        
            % Find peaks in the FFT result
            % [peaks, locs] = findpeaks(P1, 'MinPeakHeight', 0.01*max(P1));
        
            % Plot the amplitude spectrum
            figure;
            % plot(f_fft, P1);
            stem(f_fft, P1);
            hold on;
            % plot(f_fft(locs), peaks, 'ro');
            title('Single-Sided Amplitude Spectrum');
            xlabel('Frequency (Hz)');
            ylabel('Amplitude (Pa)');
            % Zoom in on the frequency range of interest
            xlim([400e3 550e3]);
            grid on;

            figure;
            stem(f_fft,phase);
            title("One-sided Phase Spectrum of Focus Signal 1");
            xlabel('Frequency(Hz)');
            ylabel('Phase (rad)')
             % Zoom in on the frequency range of interest
            xlim([400e3 550e3]);
            grid on;                             
        end
        end

% ----------------------------------------------------------------
% >>>>> Option 2: Constructing signals using createCWSignals <<<<<
% ----------------------------------------------------------------
    
    



% --------------------------------------------------
% Test effects of different params on pulse duration
% --------------------------------------------------

if plot_pulse_durations
    % f_test = [230:10:270] * 1e3;
    % T_m_fmp = find_mixed_period(f_test)

    test_ctr_freq = false;
    test_bandwidth_size = false;
    test_freq_spacing = true;

    if test_freq_spacing
        ctr_freq = 500; % kHz
        nfreq = 5; % Must be odd. Changing this paramater does not affect the graph (as long as it is odd)
        freq_spacing = [1:50];
        T_m_array = zeros(1,length(freq_spacing));
        for i = 1:length(freq_spacing)        
            lowest_freq = ctr_freq - (((nfreq - 1) / 2) * freq_spacing(i));
            highest_freq = ctr_freq + (((nfreq - 1) / 2) * freq_spacing(i));
            f0 = [lowest_freq:freq_spacing(i):highest_freq]*1e3;
            T_m = find_mixed_period(f0);
            T_m_array(i) = T_m;
        end
        figure;
        plot(freq_spacing,T_m_array * 1e3);
        title("Pulse Duration vs Frequency Spacing");
        ylabel("Pulse Duration (ms)");
        xlabel("Frequency Spacing (kHz)");          
    end

    if test_ctr_freq
        ctr_freq_array = [200:1:500];
        T_m_array = zeros(1,length(ctr_freq_array));
        for i = 1:length(ctr_freq_array)
            ctr_freq = ctr_freq_array(i); % kHz
            lowest_freq = ctr_freq - 20;
            highest_freq = ctr_freq + 20;
            f0 = [lowest_freq:10:highest_freq] * 1e3; % Hz
            T_m = find_mixed_period(f0);
            T_m_array(i) = T_m;
        end

        plot(ctr_freq_array,T_m_array);
        title("Pulse Duration based on Center Frequency");
        ylabel("Pulse Duration (s)");
        xlabel("Frequency (kHz)");
        % ylim([0 0.1])
    end

    if test_bandwidth_size
        ctr_freq = 470; % kHz    
        nfreq = [1:2:100]; % should be odd
        freq_spacing = 5;
        T_m_array = zeros(1,length(nfreq));
        for i = 1:length(nfreq)
            if i == 1
                f0 = ctr_freq * 1e3;
                T_m = find_mixed_period(f0);
            else
                lowest_freq = ctr_freq - (((nfreq(i) - 1) / 2) * freq_spacing);
                highest_freq = ctr_freq + (((nfreq(i) - 1) / 2) * freq_spacing);
                f0 = [lowest_freq:freq_spacing:highest_freq]*1e3;
                T_m = find_mixed_period(f0);
            end
            T_m_array(i) = T_m;
        end
        plot(nfreq,T_m_array);
        title("Pulse Duration based on Bandwidth Size");
        ylabel("Pulse Duration (s)");
        xlabel("Bandwidth size");
        % ylim([0 0.1])
    end
end
end




% ----------------
% Helper functions
% ----------------

function T_m = find_mixed_period(f0)
    T_array = 1.0 ./ f0;
    T_max = max(T_array);
    n = 1.0;
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