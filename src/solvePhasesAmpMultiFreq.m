function p = solvePhasesAmpMultiFreq(A, b, domain_ids, skull_ids, vol_ids, p_init, init_ids, beta)
p_init = double(p_init);

% disp('size of p_init')
% size(p_init)

% Obtain excitation vector p to have optimal amplitudes and phases to
% minimize cost function.
% 
% % Args:
%   A - cell array containing propagation matrix for each frequency
%   b - vector with desired pressures
%   p_init - initial excitation vector solution for each frequency
%
% Outputs:
%   p - optimal excitation vector

% Old code:
% Separate A and b
% [A1, A2, b1, b2, ~, ~] = prepare_opt_vars(A, b, domain_ids | skull_ids, vol_ids, init_ids);
% clear A;
% 
% % Add regularization
% [A1, b1] = add_L2_reg(A1, b1, beta(1));
% [A2, b2] = add_ineq(A2, b2, length(p_init));


nfreq = size(A, 2);

% Prepare opt vars, separate A and b into target region (A1,b1) and other indices (A2,b2)
A1_cells = {};
A2_cells = {};
for i = 1:nfreq
    A_temp = A{i};
    [A1, A2] = prepare_opt_A(A_temp, domain_ids | skull_ids, vol_ids, init_ids);
    % A1 = double(A_temp(domain_ids | skull_ids &  init_ids, :)); % Init points
    % A2 = double(A_temp(domain_ids | skull_ids & ~init_ids, :)); % Non-init points
    A1_cells{end + 1} = A1;
    A2_cells{end + 1} = A2;
end
[b1, b2] = prepare_opt_b(b, domain_ids | skull_ids, vol_ids, init_ids);


% Add regularization
for i = 1:nfreq
    [A1_cells{i}, b1] = add_L2_reg(A1_cells{i}, b1, beta(1));
    [A2_cells{i}, b2] = add_ineq(A2_cells{i}, b2, length(p_init));
end


% Cost fctn and constraints
fun = @(p)cost_fctn(p, A1_cells, b1, nfreq);
nonlcon = @(p)ineq_const(p, A2_cells, b2, nfreq);
term_fctn = @(x, optimValues, state)customOutputFcn(x, optimValues, state, 1e0, 1e0);


options = optimoptions('fmincon','Display','iter', 'FunctionTolerance', 1e-6, 'ConstraintTolerance', 1e-6, ...
    'Algorithm','active-set', 'OutputFcn', term_fctn); % interior-point, sqp, trust-region-reflective, active-set
options.MaxFunctionEvaluations = 2.5e5;
options.MaxIterations = 1e3;

[p_opt, fval, exitflag, output] = fmincon(fun, [real(p_init); imag(p_init)], [], [], [], [], [], [], nonlcon, options);

p = getCompVec(p_opt); % What does the getCompVec function do?
% The matlab optimizers cannot deal with complex numbers as input variables
% which is why we always need to separate the real and imaginary part when
% we pass p (see fmincon in line 22 -> the init vector is separated in real
% and imaginary part and then stacked). That's why every time we dive into
% the cost function and after the optimization, we need to reconvert p into
% a complex number.

end

function stop = customOutputFcn(x, optimValues, state, fval_tol, constr_tol)

    % disp("entered term func")

    stop = false;
    % Check if the absolute function value is less than the threshold
    if abs(optimValues.fval) < fval_tol && optimValues.constrviolation < constr_tol
        stop = true;
        disp('Terminating: Absolute function value is below the threshold.');
    end
end


function val = cost_fctn(p, A1_cells, b1, nfreq)

% disp("entered cost func")
% 
% disp('size of p pre getcompvec')
% disp(size(p))

p = getCompVec(p);

% disp('size of p post getcompvec')
% disp(size(p))

% Transform resulting signals into time domain and sum signals
y1 = timeDomainSum(nfreq,A1_cells,p);
% disp('size of y1')
% disp(size(y1))

% Calc cost fctn value
val = norm(abs(y1) - b1);

% disp("obtain cost func val:")
% disp(val)

end

function [c,ceq] = ineq_const(p, A2_cells, b2, nfreq)

% disp("enter ineq const")

p = getCompVec(p);

% Transform resulting signals into time domain and sum signals
y2 = timeDomainSum(nfreq,A2_cells,p);

% Calc inequality constraint value (if c > 0, then constraint is not satisfied)
c = abs(y2) - b2;
ceq = [];

% disp("obtained number of constraint violations:")
% length(find(c>0)) % returns 0 f

end



function p = getCompVec(p)

% Need to adjust this to account for multifreq (multiple columns)
p = p(1:round(end/2), :) + 1j * p(round(end/2)+1:end, :);

end

% function y = timeDomainSum(nfreq,A_cells,p)
%     n = size(A_cells{1},1);
%     y = zeros(n,1);
%     for i = 1:nfreq
%         y_temp = ifft(A_cells{i} * p(:,i));
%         y = [y + y_temp];
%     end
% end
