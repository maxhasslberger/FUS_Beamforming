function y = timeDomainSum(f0,A_cells,p)
    n = size(A_cells{1},1);
    nfreq = length(f0);
    T = 1/(min(f0));
    Fs = 100e6;
    dt = 1/Fs;
    t = 0:dt:10*T;
    signals = zeros(n,length(t));
    for i = 1:nfreq
        w = 2 * pi * f0(i);
        y_temp = A_cells{i} * p(:,i);
        amp = abs(y_temp); phi = angle(y_temp);
        % signals_tmp = amp .* cos(w .* t + phi);        
        % signals = [signals + signals_tmp];
        for j = 1:n
            signals(j, :) = signals(j, :) + amp(j) * cos(w * t + phi(j));
        end
    end
    y = (max(signals, [], 2) - min(signals, [], 2)) ./ 2;


    % % Old code:
    % n = size(A_cells{1},1);
    % y = zeros(n,1);
    % nfreq = length(f0);
    % for i = 1:nfreq
    %     y_temp = ifft(A_cells{i} * p(:,i));
    %     y = [y + y_temp];
    % end
end
