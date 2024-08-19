function p = solvePhasesAmp(A, b, cons_ids, vol_ids, p_init, init_ids, beta, ineq_active)
p_init = double(p_init);

% Separate A and b
[A1, A2, b1, b2] = prepare_opt_vars(A, b, cons_ids, vol_ids, init_ids, ineq_active);

clear A;

A0 = zeros(1, length(p_init));
b0 = 0.0;

% Add regularization
[A0, b0] = add_L2_reg(A0, b0, beta(1));

% Optimization options
term_fctn = @(x, optimValues, state)customOutputFcn(x, optimValues, state, [], 1e-3);
options = optimoptions('fmincon','Display','iter', 'FunctionTolerance', 1e6, 'ConstraintTolerance', 1e-3, ...
    'Algorithm','active-set');%, 'OutputFcn', term_fctn); % interior-point, sqp, trust-region-reflective, active-set
options.MaxFunctionEvaluations = 2.5e5;
options.MaxIterations = 1e3;

% Define cost function
fun = @(p)cost_fctn(p, A0, b0);

% Optimize
if false
    [p, fval, exitflag, output] = solve_ip(fun, A1, b1, A2, b2, p_init, options);
else
    A2_iter = zeros(1, length(p_init));
    b2_iter = 0.0;
    p = p_init;

    max_iter = 20;
    tol = b2 / 1000; % Pa
    for iter = 1:max_iter
        disp(strcat("Main iteration ", num2str(iter), "..."))

        % Optimize and check if ineq constraints fulfilled
        [p, fval, exitflag, output] = solve_ip(fun, A1, b1, A2_iter, b2_iter, p, options);
        b2_new = abs(A2 * p);

        unful_ineq = b2_new > (b2 + tol);
        if ~any(unful_ineq)
            break; % successful
        else
            % Add inequalities that are not fulfilled
            A2_iter = [A2_iter; A2(unful_ineq, :)];
            b2_iter = [b2_iter; b2(unful_ineq)];
            disp(strcat("Violated Inequalities: ", num2str(sum(unful_ineq))))
        end
    end

    if any(unful_ineq)
        disp("Optimization failed");
    else
        disp("Optimization successful");
    end
end

end

function [p, fval, exitflag, output] = solve_ip(fun, A1, b1, A2, b2, p_init, options)

% Update constraints
nonlcon = @(p)ineq_const(p, A1, b1, A2, b2);

% Optimize
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
