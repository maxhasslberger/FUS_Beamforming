function p = solvePhasesOnly(A, b, cons_ids, vol_ids, p_init, init_ids, beta, ineq_active, mask2el, el_per_t, via_abs)
p_init = double(p_init);

% Separate A and b
[A1, A2, b1, b2, ~, ~] = prepare_opt_vars(A, b, cons_ids, vol_ids, init_ids, ineq_active);

clear A;

% Add regularization and constraints
[A2, b2] = add_ineq(A2, b2, length(p_init));

% Obtain sub matrices
A1_r = real(A1);
A1_i = imag(A1);
A2_r = real(A2);
A2_i = imag(A2);

% Define one amp per transducer
n_amps = length(el_per_t);
trx_ids = cell(1, n_amps);
shift = 0;
for i = 1:n_amps
    trx_ids{i} = mask2el(1 + shift:sum(el_per_t(1:i)));
    shift = shift + el_per_t(i);
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
    fun = @(p)cost_fctn(p, A1_r, A1_i, b1, trx_ids);
    nonlcon = @(p)unitdisk(p, A1_r, A1_i, b1, A2_r, A2_i, b2, via_abs, amp_fac, trx_ids);
else
    % Define initial vector
    p_start = zeros(length(p_init) + n_amps, 1);
    amp_fac = init_amp;
    p_start(1:n_amps) = 1.0;

    p_start(n_amps + 1:end) = angle(p_init);

    % Cost fctn and constraints
    fun = @(p)cost_fctn2(p, A1, b1, amp_fac, trx_ids);
    nonlcon = @(p)unitdisk(p, A1_r, A1_i, b1, A2_r, A2_i, b2, via_abs, amp_fac, trx_ids);
end
term_fctn = @(x, optimValues, state)customOutputFcn(x, optimValues, state, 1e0, 1e0);

options = optimoptions('fmincon','Display','iter', 'FunctionTolerance', 1e-6, 'ConstraintTolerance', 1e-6, ...
    'Algorithm','active-set', 'OutputFcn', term_fctn); % interior-point, sqp, trust-region-reflective, active-set
options.MaxFunctionEvaluations = 2.5e5;
options.MaxIterations = 1e3;


[p_opt, fval, exitflag, output] = fmincon(fun, p_start, [], [], [], [], [], [], nonlcon, options);


p_0 = zeros(size(A1, 2), 1);
p_0 = getAmpPerElement(p_0, p_opt, trx_ids);

if via_abs
    p = p_0 .* exp(1j * atan2(p_opt(n_amps + ceil((end-n_amps)/2) + 1:end), p_opt(n_amps + (1:ceil((end-n_amps)/2)))));
else
    p = amp_fac * p_0 .* exp(1j * p_opt(n_amps + 1:end));
end


end

function stop = customOutputFcn(x, optimValues, state, fval_tol, constr_tol)
    stop = false;
    % Check if the absolute function value is less than the threshold
    if abs(optimValues.fval) < fval_tol && optimValues.constrviolation < constr_tol
        stop = true;
        disp('Terminating: Absolute function value is below the threshold.');
    end
end


function val = cost_fctn(p, A1_r, A1_i, b1, trx_ids)

% [x_r, x_i] = getElements_abs(p, trx_ids);
% val = norm(sqrt((A1_r * x_r - A1_i * x_i).^2 + (A1_i * x_r + A1_r * x_i).^2) - b1);
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
