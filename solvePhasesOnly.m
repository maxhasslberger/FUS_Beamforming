
function p = solvePhasesOnly(A1, A2, b1, b2, via_abs)

A1 = double(A1);
A2 = double(A2);
% Obtain sub matrices
A1_r = real(A1);
A1_i = imag(A1);
if ~isempty(A2)
    A2_r = real(A2);
    A2_i = imag(A2);
else
    A2_r = zeros(1, size(A1, 2));
    A2_i = zeros(1, size(A1, 2));
end

% Cost fctn and constraints
fun = @(p)cost_fctn(p, A1_r, A1_i, b1, via_abs);
nonlcon = @(p)unitdisk(p, A2_r, A2_i, b2, via_abs);

options = optimoptions('fmincon','Display','iter', 'FunctionTolerance', 1e-6, 'ConstraintTolerance', 1e-6, ...
    'Algorithm','active-set'); % interior-point, sqp, trust-region-reflective, active-set

% Define initial solution
p_lin = pinv(A1) * b1;

if via_abs
    p_start = zeros(length(p_lin) * 2 + 1, 1);
    p_start(1) = mean(abs(p_lin));

    p_start(2:ceil(end/2)) = real(p_lin);
    p_start(ceil(end/2) + 1:end) = imag(p_lin); % TODO: Numerical error bc of div by 0?
else
    p_start = zeros(length(p_lin) + 1, 1);
    p_start(1) = mean(abs(p_lin));

    p_start(2:end) = angle(p_lin); % TODO: Numerical precision bc of first element...
end

[p_opt, fval, exitflag, output] = fmincon(fun, p_start, [], [], [], [], [], [], nonlcon, options);


p = p_opt(1) * exp(1j * atan2(p_opt(ceil(end/2) + 1:end), p_opt(2:ceil(end/2))));

end


function val = cost_fctn(p, A1_r, A1_i, b1, via_abs)

if via_abs
    [x_r, x_i, p_0] = getElements_abs(p);
else
    [x_r, x_i, p_0] = getElements_phase(p);
end


val = sum((p_0 * sqrt((A1_r * x_r - A1_i * x_i).^2 + (A1_i * x_r + A1_r * x_i).^2) - b1).^2);

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
