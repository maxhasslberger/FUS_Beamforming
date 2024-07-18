function y = timeDomainSum(nfreq,A_cells,p)
    n = size(A_cells{1},1);
    y = zeros(n,1);
    for i = 1:nfreq
        y_temp = ifft(A_cells{i} * p(:,i));
        y = [y + y_temp];
    end
end
