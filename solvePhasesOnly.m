function p = solvePhasesOnly(A, b, opt_ids, obs_ids, p_init, init_ids, beta, mask2el, el_per_t, via_abs)
p_init = double(p_init);

% Separate A and b
[A1, A2, b1, b2, A_zero, A_vol, gamma] = prepare_opt_vars(A, b, opt_ids, obs_ids, init_ids);
clear A;

% Add regularization
[A1, b1, gamma] = add_L2_reg(A1, b1, gamma, beta(1));
[A1, b1, gamma] = add_zeroAmp_reg(A1, b1, gamma, A_zero, beta(2));
[A1, b1, gamma] = add_volAmp_reg(A1, b1, gamma, A_vol, beta(3));
[A2, b2] = add_ineq(A2, b2, length(p_init), beta(4));

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
    fun = @(p)cost_fctn(p, A1_r, A1_i, b1, gamma, trx_ids);
    nonlcon = @(p)unitdisk(p, A2_r, A2_i, b2, via_abs, amp_fac, trx_ids);
else
    % Define initial vector
    p_start = zeros(length(p_init) + n_amps, 1);
    amp_fac = init_amp;
    p_start(1:n_amps) = 1.0;

    p_start(n_amps + 1:end) = angle(p_init);

    % Cost fctn and constraints
    fun = @(p)cost_fctn2(p, A1, b1, gamma, amp_fac, trx_ids);
    nonlcon = @(p)unitdisk(p, A2_r, A2_i, b2, via_abs, amp_fac, trx_ids);
end

options = optimoptions('fmincon','Display','iter', 'FunctionTolerance', 1e-6, 'ConstraintTolerance', 1e-6, ...
    'Algorithm','active-set'); % interior-point, sqp, trust-region-reflective, active-set
options.MaxFunctionEvaluations = 0.1e5;


[p_opt, fval, exitflag, output] = fmincon(fun, p_start, [], [], [], [], [], [], nonlcon, options);


p_0 = zeros(size(A1, 2), 1);
p_0 = getAmpPerElement(p_0, p_opt, trx_ids);

if via_abs
    p = p_0 .* exp(1j * atan2(p_opt(n_amps + ceil((end-n_amps)/2) + 1:end), p_opt(n_amps + (1:ceil((end-n_amps)/2)))));
else
    p = amp_fac * p_0 .* exp(1j * p_opt(n_amps + 1:end));
end


end


function val = cost_fctn(p, A1_r, A1_i, b1, gamma, trx_ids)

[x_r, x_i] = getElements_abs(p, trx_ids);
val = norm(sqrt((A1_r * x_r - A1_i * x_i).^2 + (A1_i * x_r + A1_r * x_i).^2).^gamma - b1);

end

function val = cost_fctn2(p, A1, b1, gamma, amp_fac, trx_ids)
n_amps = length(trx_ids);

p_0 = zeros(size(A1, 2), 1);
p_0 = getAmpPerElement(p_0, p, trx_ids);

val = norm(amp_fac * abs(A1 * (p_0 .* exp(1j * p(n_amps + 1:end)))).^gamma - b1);

end


function [x_r, x_i] = getElements_abs(p, trx_ids)

n_amps = length(trx_ids);

p_r = p(n_amps + (1:ceil((end-n_amps)/2)));
p_i = p(n_amps + ceil((end-n_amps)/2) + 1:end);

p_0 = zeros(size(p_r));
p_0 = getAmpPerElement(p_0, p, trx_ids);

x_r = p_0 .* p_r ./ sqrt(p_r.^2 + p_i.^2);
x_i = p_0 .* p_i ./ sqrt(p_r.^2 + p_i.^2);

end

function p_0 = getAmpPerElement(p_0, p, trx_ids)

n_amps = length(trx_ids);
for i = 1:n_amps
    p_0(trx_ids{i}) = p(i);
end

end
