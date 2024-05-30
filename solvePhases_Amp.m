function p = solvePhases_Amp(A, b, opt_ids, obs_ids, p_init, init_ids, beta)
p_init = double(p_init);

% Separate A and b
[A1, A2, b1, b2, A_zero, A_vol, gamma] = prepare_opt_vars(A, b, opt_ids, obs_ids, init_ids);
clear A;

% Add regularization
[A1, b1, gamma] = add_L2_reg(A1, b1, gamma, beta(1));
[A1, b1, gamma] = add_zeroAmp_reg(A1, b1, gamma, A_zero, beta(2));
[A1, b1, gamma] = add_volAmp_reg(A1, b1, gamma, A_vol, beta(3));
[A2, b2] = add_ineq(A2, b2, length(p_init), beta(4));

% Cost fctn and constraints
fun = @(p)cost_fctn(p, A1, b1, gamma);
nonlcon = @(p)ineq_const(p, A2, b2);
term_fctn = @(x, optimValues, state)customOutputFcn(x, optimValues, state, 1e0, 1e0);

options = optimoptions('fmincon','Display','iter', 'FunctionTolerance', 1e-6, 'ConstraintTolerance', 1e-6, ...
    'Algorithm','active-set', 'OutputFcn', term_fctn); % interior-point, sqp, trust-region-reflective, active-set
options.MaxFunctionEvaluations = 0.5e5;

[p_opt, fval, exitflag, output] = fmincon(fun, [real(p_init); imag(p_init)], [], [], [], [], [], [], nonlcon, options);

p = getCompVec(p_opt);

end

function stop = customOutputFcn(x, optimValues, state, fval_tol, constr_tol)
    stop = false;
    % Check if the absolute function value is less than the threshold
    if abs(optimValues.fval) < fval_tol && optimValues.constrviolation < constr_tol
        stop = true;
        disp('Terminating: Absolute function value is below the threshold.');
    end
end


function val = cost_fctn(p, A1, b1, gamma)

p = getCompVec(p);
val = norm(abs(A1 * p).^gamma - b1);

end

function [c,ceq] = ineq_const(p, A2, b2)

p = getCompVec(p);
c = abs(A2 * p) - b2;

ceq = [];

end

function p = getCompVec(p)

p = p(1:round(end/2)) + 1j * p(round(end/2)+1:end);

end
