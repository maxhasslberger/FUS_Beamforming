function p = solvePhases_Amp(A, b, p_init, init_ids, beta_L2)

% Separate A and b
A1 = double(A(init_ids, :));
A2 = double(A);
clear A;
A2(init_ids, :) = [];

b1 = double(b(init_ids));
b2 = double(b);
clear b;
b2(init_ids) = [];

p_init = double(p_init);

% Add L2 regularization
[A1, b1] = add_L2_reg(A1, b1, beta_L2);

if isempty(A2)
    A2 = zeros(1, size(A1, 2));
end

% Cost fctn and constraints
fun = @(p)cost_fctn(p, A1, b1);
nonlcon = @(p)ineq_const(p, A2, b2);

options = optimoptions('fmincon','Display','iter', 'FunctionTolerance', 1e-6, 'ConstraintTolerance', 1e-6, ...
    'Algorithm','active-set'); % interior-point, sqp, trust-region-reflective, active-set
options.MaxFunctionEvaluations = 0.5e5;

[p_opt, fval, exitflag, output] = fmincon(fun, [real(p_init); imag(p_init)], [], [], [], [], [], [], nonlcon, options);

p = getCompVec(p_opt);

end


function val = cost_fctn(p, A1, b1)

p = getCompVec(p);
val = sum((abs(A1 * p) - b1).^2);

end

function [c,ceq] = ineq_const(p, A2, b2)

p = getCompVec(p);
c = abs(A2 * p) - b2;

ceq = [];

end

function p = getCompVec(p)

p = p(1:round(end/2)) + 1j * p(round(end/2)+1:end);

end
