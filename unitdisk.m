function [c,ceq] = unitdisk(p, A2_r, A2_i, b2, via_abs, amp_fac)

if via_abs
    [x_r, x_i, p_0] = getElements_abs(p);
    c = (p_0 * sqrt((A2_r * x_r - A2_i * x_i).^2 + (A2_i * x_r + A2_r * x_i).^2) - b2);
else
    c = sum((p(1) * amp_fac * abs((A2_r + 1j * A2_i) * exp(1j*p(2:end))) - b2));
end

ceq = [];

end


function [x_r, x_i, p_0] = getElements_abs(p)

p_r = p(2:ceil(end/2));
p_i = p(ceil(end/2) + 1:end);
p_0 = p(1);

x_r = p_r ./ sqrt(p_r.^2 + p_i.^2);
x_i = p_i ./ sqrt(p_r.^2 + p_i.^2);

end
