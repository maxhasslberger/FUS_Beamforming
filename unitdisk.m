function [c,ceq] = unitdisk(p, A2_r, A2_i, b2, via_abs)

if via_abs
    [x_r, x_i, p_0] = getElements_abs(p);
else
    [x_r, x_i, p_0] = getElements_phase(p);
end

c = (p_0 * sqrt((A2_r * x_r - A2_i * x_i).^2 + (A2_i * x_r + A2_r * x_i)) - b2);
ceq = [];

end


function [x_r, x_i, p_0] = getElements_abs(p)

p_r = p(2:ceil(end/2));
p_i = p(ceil(end/2) + 1:end);
p_0 = p(1);

x_r = p_r ./ sqrt(p_r.^2 + p_i.^2);
x_i = p_i ./ sqrt(p_r.^2 + p_i.^2);

end

function [x_r, x_i, p_0] = getElements_phase(p)

phase_vec = exp(1j * p(2:end));
p_0 = p(1);

x_r = real(phase_vec);
x_i = imag(phase_vec);

end
