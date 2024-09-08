function p = solvePhasesAmp(iter_mode, A1, A2, b1, b2, p_init, beta, term)
p_init = double(p_init);


A0 = zeros(1, length(p_init));
b0 = 0.0;

% Add regularization
[A0, b0] = add_L2_reg(A0, b0, beta(1));

if isempty(term)
    term.fun_tol = 1e-1;
    term.constr_tol = 1e-3;
    term.iter_tol = 10;
    term.iter_lim = max([term.iter_tol, 200 - term.iter_tol]);
    term.norm_val = 10;
end

% Optimization options
term_fctn = @(x, optimValues, state)customOutputFcn(x, optimValues, state, term.fun_tol, term.constr_tol, term.iter_lim);
options = optimoptions('fmincon','Display','iter', 'FunctionTolerance', term.fun_tol, 'ConstraintTolerance', term.constr_tol, ...
    'Algorithm', 'active-set', 'MaxIterations', term.iter_lim + term.iter_tol, 'MaxFunctionEvaluations', inf, 'OutputFcn', term_fctn); 
% Algorithms: interior-point, sqp, trust-region-reflective, active-set

% Define cost function
fun = @(p)cost_fctn(p, A0, b0, term.norm_val);

% Optimize
if ~iter_mode
    [p, fval, exitflag, output] = solve_ip(fun, A1, b1, A2, b2, p_init, options);
    if exitflag == -2 % No feasible point found
        disp(output.message)
        disp("Optimization failed");
    else
        disp("Optimization successful");
    end
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
        if exitflag == -2 % No feasible point found
            unful_ineq = 1;
            disp(output.message)
            break;
        end

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

function stop = customOutputFcn(x, optimValues, state, fun_tol, constr_tol, iter_lim)

stop = false;
persistent prevFval % Store the previous function value
currFval = optimValues.fval;

if strcmp(state, 'init')
    % Initialize prevFval during the first iteration
    prevFval = inf;
end
        
% Check if the values are less than the thresholds
fun_tol_curr = (prevFval - currFval) / (currFval + 1);
if optimValues.constrviolation < constr_tol && (currFval == 0 || (fun_tol_curr < fun_tol && fun_tol_curr > 0) || optimValues.iteration >= iter_lim)
    stop = true;
    disp(strcat("Terminating: Optimality conditions fulfilled. Rel. function tolerance: ", ...
        num2str(fun_tol_curr), ", Max. constraint violation: ", num2str(optimValues.constrviolation)));
end

% Update prevFval for the next iteration
prevFval = currFval;

end


function val = cost_fctn(p, A0, b0, norm_val)

p = getCompVec(p);
val = norm(abs(A0 * p) - b0, norm_val);

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
