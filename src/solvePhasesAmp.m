function p = solvePhasesAmp(A, b, cons_ids, vol_ids, p_init, init_ids, beta, ineq_active)
p_init = double(p_init);

if ineq_active
    % Separate A and b
    [A1, A2, b1, b2, ~, ~] = prepare_opt_vars(A, b, cons_ids, vol_ids, init_ids);
else
    A1 = double(A(init_ids, :));
    b1 = double(b(init_ids));
    A2 = [];
    b2 = [];
end

clear A;

A0 = zeros(1, length(p_init));
b0 = 0.0;

% Add regularization
[A0, b0] = add_L2_reg(A0, b0, beta(1));
[A2, b2] = add_ineq(A2, b2, length(p_init));

% Cost fctn and constraints
fun = @(p)cost_fctn(p, A0, b0);
nonlcon = @(p)ineq_const(p, A1, b1, A2, b2);
term_fctn = @(x, optimValues, state)customOutputFcn(x, optimValues, state, [], 1e-3);

options = optimoptions('fmincon','Display','iter', 'FunctionTolerance', 1e6, 'ConstraintTolerance', 1e-3, ...
    'Algorithm','active-set');%, 'OutputFcn', term_fctn); % interior-point, sqp, trust-region-reflective, active-set
options.MaxFunctionEvaluations = 2.5e5;
options.MaxIterations = 1e3;

[p_opt, fval, exitflag, output] = fmincon(fun, [real(p_init); imag(p_init)], [], [], [], [], [], [], nonlcon, options);

p = getCompVec(p_opt);

end

function stop = customOutputFcn(x, optimValues, state, fval_tol, constr_tol)
    stop = false;
    % Check if the absolute function value is less than the threshold
    if optimValues.constrviolation < constr_tol % && abs(optimValues.fval) < fval_tol
        stop = true;
        disp('Terminating: Absolute function value is below the threshold.');
    end
end


function val = cost_fctn(p, A0, b0)

p = getCompVec(p);
val = norm(abs(A0 * p) - b0);

end

function [c,ceq] = ineq_const(p, A1, b1, A2, b2)

p = getCompVec(p);

% ceq = [];
ceq = abs(A1 * p) - b1;
c = abs(A2 * p) - b2;

end

function p = getCompVec(p)

p = p(1:round(end/2)) + 1j * p(round(end/2)+1:end);

end
