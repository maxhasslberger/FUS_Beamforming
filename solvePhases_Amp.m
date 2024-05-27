function p = solvePhases_Amp(A, b, opt_ids, obs_ids, p_init, init_ids, beta)

% Separate A and b
[A1, A2, b1, b2, b_zero, b_vol, p_init] = prepare_opt_vars(A, b, p_init, opt_ids, obs_ids, init_ids);
clear A;

% Add L2 regularization
[A1, b1] = add_L2_reg(A1, b1, beta(1));
[A1, b1] = add_zeroAmp_reg(A1, b1, b_zero, beta(2));
[A1, b1] = add_volAmp_reg(A1, b1, b_vol, beta(3));

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
