function [c,ceq] = unitdisk(p, A2_r, A2_i, b2, via_abs, amp_fac, trx_ids)

if via_abs
    [x_r, x_i] = getElements_abs(p, trx_ids);
    c = sqrt((A2_r * x_r - A2_i * x_i).^2 + (A2_i * x_r + A2_r * x_i).^2) - b2;
else
    n_amps = length(trx_ids);

    p_0 = zeros(size(A2, 2), 1);
    p_0 = getAmpPerElement(p_0, p, trx_ids);

    c = amp_fac * abs((A2_r + 1j * A2_i) * (p_0 .* exp(1j * p(n_amps:end)))) - b2;
end

ceq = [];

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
