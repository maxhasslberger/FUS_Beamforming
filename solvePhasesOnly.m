
function p = solvePhasesOnly(A1, A2, b1, b2, p_init, beta_L2, mask2el, el_per_t, via_abs)

A1 = double(A1);
A2 = double(A2);

% Add L2 regularization
[A1, b1] = add_L2_reg(A1, b1, beta_L2);

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

% Define one amp per transducer
n_amps = length(el_per_t);
trx_ids = cell(1, n_amps);
shift = 0;
for i = 1:n_amps
    trx_ids{i} = mask2el(1 + shift:sum(el_per_t(1:i)));
    shift = shift + el_per_t(i);
end

p_init = double(p_init);
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
    nonlcon = @(p)unitdisk(p, A2_r, A2_i, b2, via_abs, amp_fac, trx_ids);
else
    % Define initial vector
    p_start = zeros(length(p_init) + n_amps, 1);
    amp_fac = init_amp;
    p_start(1:n_amps) = 1.0;

    p_start(n_amps + 1:end) = angle(p_init);

    % Cost fctn and constraints
    fun = @(p)cost_fctn2(p, A1, b1, amp_fac, trx_ids);
    nonlcon = @(p)unitdisk(p, A2_r, A2_i, b2, via_abs, amp_fac, trx_ids);
end

options = optimoptions('fmincon','Display','iter', 'FunctionTolerance', 1e-6, 'ConstraintTolerance', 1e-6, ...
    'Algorithm','active-set'); % interior-point, sqp, trust-region-reflective, active-set
options.MaxFunctionEvaluations = 0.5e5;


[p_opt, fval, exitflag, output] = fmincon(fun, p_start, [], [], [], [], [], [], nonlcon, options);


p_0 = zeros(size(A1, 2), 1);
p_0 = getAmpPerElement(p_0, p_opt, trx_ids);

if via_abs
    p = p_0 .* exp(1j * atan2(p_opt(n_amps + ceil((end-n_amps)/2) + 1:end), p_opt(n_amps + (1:ceil((end-n_amps)/2)))));
else
    p = amp_fac * p_0 .* exp(1j * p_opt(n_amps + 1:end));
end


end


function val = cost_fctn(p, A1_r, A1_i, b1, trx_ids)

[x_r, x_i] = getElements_abs(p, trx_ids);
val = sum((sqrt((A1_r * x_r - A1_i * x_i).^2 + (A1_i * x_r + A1_r * x_i).^2) - b1).^2);

end

function val = cost_fctn2(p, A1, b1, amp_fac, trx_ids)
n_amps = length(trx_ids);

p_0 = zeros(size(A1, 2), 1);
p_0 = getAmpPerElement(p_0, p, trx_ids);

val = sum((amp_fac * abs(A1 * (p_0 .* exp(1j * p(n_amps + 1:end)))) - b1).^2);

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
