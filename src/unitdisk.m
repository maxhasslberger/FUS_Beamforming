function [c,ceq] = unitdisk(p, A1, b1, A2, b2, via_abs, amp_fac, trx_ids)

if via_abs
    [x_r, x_i] = getElements_abs(p, trx_ids);
    c = sqrt((real(A2) * x_r - imag(A2) * x_i).^2 + (imag(A2) * x_r + real(A2) * x_i).^2) - b2;
    ceq = sqrt((real(A1) * x_r - imag(A1) * x_i).^2 + (imag(A1) * x_r + real(A1) * x_i).^2) - b1;
else
    n_amps = length(trx_ids);

    p_0 = zeros(size(A2, 2), 1);
    p_0 = getAmpPerElement(p_0, p, trx_ids);
    c = amp_fac * abs(A2 * (p_0 .* exp(1j * p(n_amps + 1:end)))) - b2;
    
    p_0 = zeros(size(A1, 2), 1);
    p_0 = getAmpPerElement(p_0, p, trx_ids);
    ceq = amp_fac * abs(A1 * (p_0 .* exp(1j * p(n_amps + 1:end)))) - b1;
end

end


function [x_r, x_i] = getElements_abs(p, trx_ids)

n_amps = length(trx_ids);

p_r = p(n_amps + (1:ceil((end-n_amps)/2)));
p_i = p(n_amps + ceil((end-n_amps)/2) + 1:end);

p_0 = zeros(size(p_r));
p_0 = getAmpPerElement(p_0, p, trx_ids);

x_r = p_0 .* p_r ./ sqrt(p_r.^2 + p_i.^2);
x_i = p_0 .* p_i ./ sqrt(p_r.^2 + p_i.^2);

end

function p_0 = getAmpPerElement(p_0, p, trx_ids)

n_amps = length(trx_ids);
for i = 1:n_amps
    p_0(trx_ids{i}) = p(i);
end

end
