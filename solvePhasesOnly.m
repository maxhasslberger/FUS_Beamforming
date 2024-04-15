
function p = solvePhasesOnly(A1, A2, b1, b2)

fun = @(p)cost_fctn(p, A1, b1);
nonlcon = @(p)unitdisk(p, A2, b2);
p_start = pinv(A1) * b1;

[p_opt, fval, exitflag, output] = fmincon(fun, p_start, [], [], [], [], [], [], nonlcon);


p = p_opt(1) * exp(1j * atan2(p_opt(ceil(end/2) + 1:end), p_opt(2:ceil(end/2))));

end


function val = cost_fctn(p, A1, b1)

[x_r, x_i, p_0, A1_r, A1_i] = getElements(p, A1);

val = sum((p_0 * sqrt((A1_r * x_r - A1_i * x_i).^2 + (A1_i * x_r + A1_r * x_i)) - b1).^2);

end


function [x_r, x_i, p_0, A_r, A_i] = getElements(p, A)

p_r = p(2:ceil(end/2));
p_i = p(ceil(end/2) + 1:end);
p_0 = p(1);

A_r = A(1:end/2, :);
A_i = A(end/2 + 1:end, :);

x_r = p_r ./ sqrt(p_r.^2 + p_i.^2);
x_i = p_i ./ sqrt(p_r.^2 + p_i.^2);

end
