function p = solvePhasesOnly(iter_mode, A1, A2, b1, b2, p_init, beta, term, mask2el, el_per_t, via_abs, tot_1amp)
p_init = double(p_init);

clear A;

if ~tot_1amp
    % Obtain one amp per transducer
    n_amps = length(el_per_t);
    trx_ids = cell(1, n_amps);
    shift = 0;
    for i = 1:n_amps
        trx_ids{i} = mask2el(1 + shift:sum(el_per_t(1:i)));
        shift = shift + el_per_t(i);
    end
else
    % Obtain one amp for all transducers
    n_amps = 1;
    trx_ids = {(1:sum(el_per_t))};
end

init_amp = mean(abs(p_init));

if via_abs
    % Construct initial vector
    p_start = zeros(length(p_init) * 2 + n_amps, 1);
    amp_fac = 0.0;
    p_start(1:n_amps) = init_amp;

    p_start(n_amps + (1:ceil((end-n_amps)/2))) = real(p_init);
    p_start(n_amps + ceil((end-n_amps)/2) + 1:end) = imag(p_init);

    % Cost fctn and constraints
    fun = @(p)cost_fctn(p, A1, b1, trx_ids);
    nonlcon = @(p)unitdisk(p, A1, b1, A2, b2, via_abs, amp_fac, trx_ids);
else
    % Define initial vector
    p_start = zeros(length(p_init) + n_amps, 1);
    amp_fac = init_amp;
    p_start(1:n_amps) = 1.0;

    p_start(n_amps + 1:end) = angle(p_init);

    % Cost fctn and constraints
    fun = @(p)cost_fctn2(p, A1, b1, amp_fac, trx_ids);
    nonlcon = @(p)unitdisk(p, A1, b1, A2, b2, via_abs, amp_fac, trx_ids);
end

if isempty(term)
    term.fun_tol = 1e-1;
    term.constr_tol = 1e-3;
    term.iter_tol = 10;
    term.iter_lim = max([term.iter_tol, 200 - term.iter_tol]);
    term.norm_val = 10;
end

% Optimization options
% term_fctn = @(x, optimValues, state)customOutputFcn(x, optimValues, state, term.fun_tol, term.constr_tol, term.iter_lim);
options = optimoptions('fmincon','Display','iter', 'FunctionTolerance', term.fun_tol, 'ConstraintTolerance', term.constr_tol, ...
    'Algorithm', 'active-set', 'MaxIterations', term.iter_lim + term.iter_tol, 'MaxFunctionEvaluations', inf); % , 'OutputFcn', term_fctn); 
% Algorithms: interior-point, sqp, trust-region-reflective, active-set

if ~iter_mode
    [p_opt, fval, exitflag, output] = solve_ip(fun, A1, b1, A2, b2, p_start, options, via_abs, amp_fac, trx_ids);

    p_0 = zeros(size(A1, 2), 1);
    p_0 = getAmpPerElement(p_0, p_opt, trx_ids);
    
    if via_abs
        p = p_0 .* exp(1j * atan2(p_opt(n_amps + ceil((end-n_amps)/2) + 1:end), p_opt(n_amps + (1:ceil((end-n_amps)/2)))));
    else
        p = amp_fac * p_0 .* exp(1j * p_opt(n_amps + 1:end));
    end

    if exitflag == -2 % No feasible point found
        disp(output.message)
        disp("Optimization failed");
    else
        disp("Optimization successful");
    end
else
    A2_iter = zeros(1, length(p_init));
    b2_iter = 0.0;
    p_opt = p_start;

    max_iter = 20;
    tol = b2 / 1000; % Pa
    for iter = 1:max_iter
        disp(strcat("Main iteration ", num2str(iter), "..."))

        % Optimize and check if ineq constraints fulfilled
        [p_opt, fval, exitflag, output] = solve_ip(fun, A1, b1, A2_iter, b2_iter, p_opt, options, via_abs, amp_fac, trx_ids);
        if exitflag == -2 % No feasible point found
            unful_ineq = 1;
            disp(output.message)
            break;
        end

        p_0 = zeros(size(A1, 2), 1);
        p_0 = getAmpPerElement(p_0, p_opt, trx_ids);
        p = p_0 .* exp(1j * atan2(p_opt(n_amps + ceil((end-n_amps)/2) + 1:end), p_opt(n_amps + (1:ceil((end-n_amps)/2)))));

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

% [p_opt, fval, exitflag, output] = fmincon(fun, p_start, [], [], [], [], [], [], nonlcon, options);


% p_0 = zeros(size(A1, 2), 1);
% p_0 = getAmpPerElement(p_0, p_opt, trx_ids);
% 
% if via_abs
%     p = p_0 .* exp(1j * atan2(p_opt(n_amps + ceil((end-n_amps)/2) + 1:end), p_opt(n_amps + (1:ceil((end-n_amps)/2)))));
% else
%     p = amp_fac * p_0 .* exp(1j * p_opt(n_amps + 1:end));
% end


end

function [p_opt, fval, exitflag, output] = solve_ip(fun, A1, b1, A2, b2, p_init, options, via_abs, amp_fac, trx_ids)

% Update constraints
nonlcon = @(p)unitdisk(p, A1, b1, A2, b2, via_abs, amp_fac, trx_ids);

% Optimize
[p_opt, fval, exitflag, output] = fmincon(fun, p_init, [], [], [], [], [], [], nonlcon, options);

end

function stop = customOutputFcn(x, optimValues, state, fval_tol, constr_tol)
    stop = false;
    % Check if the absolute function value is less than the threshold
    if abs(optimValues.fval) < fval_tol && optimValues.constrviolation < constr_tol
        stop = true;
        disp('Terminating: Absolute function value is below the threshold.');
    end
end


function val = cost_fctn(p, A1, b1, trx_ids)

% [x_r, x_i] = getElements_abs(p, trx_ids);
% val = norm(sqrt((real(A1) * x_r - imag(A1) * x_i).^2 + (imag(A1) * x_r + real(A1) * x_i).^2) - b1);
val = 0;

end

function val = cost_fctn2(p, A1, b1, amp_fac, trx_ids)
% n_amps = length(trx_ids);
% 
% p_0 = zeros(size(A1, 2), 1);
% p_0 = getAmpPerElement(p_0, p, trx_ids);
% 
% val = norm(amp_fac * abs(A1 * (p_0 .* exp(1j * p(n_amps + 1:end)))) - b1);
val = 0;

end


% function [x_r, x_i] = getElements_abs(p, trx_ids)
% 
% n_amps = length(trx_ids);
% 
% p_r = p(n_amps + (1:ceil((end-n_amps)/2)));
% p_i = p(n_amps + ceil((end-n_amps)/2) + 1:end);
% 
% p_0 = zeros(size(p_r));
% p_0 = getAmpPerElement(p_0, p, trx_ids);
% 
% x_r = p_0 .* p_r ./ sqrt(p_r.^2 + p_i.^2);
% x_i = p_0 .* p_i ./ sqrt(p_r.^2 + p_i.^2);
% 
% end

function p_0 = getAmpPerElement(p_0, p, trx_ids)

n_amps = length(trx_ids);
for i = 1:n_amps
    p_0(trx_ids{i}) = p(i);
end

end
