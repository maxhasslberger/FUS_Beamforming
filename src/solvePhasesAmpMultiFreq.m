function p = solvePhasesAmpMultiFreq(A, b, domain_ids, skull_ids, vol_ids, p_init, init_ids, beta, f0)
p_init = double(p_init);

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

A0 = zeros(1, size(p_init, 1));
b0 = zeros(1, size(p_init, 2));


% Add regularization
for i = 1:nfreq
    [A1_cells{i}, b1] = add_L2_reg(A1_cells{i}, b1, beta(1));
    [A2_cells{i}, b2] = add_ineq(A2_cells{i}, b2, length(p_init));
end
[A0, b0] = add_L2_reg(A0, b0, beta(1));

term_fctn = @(x, optimValues, state)customOutputFcn(x, optimValues, state, 1e0, 1e0);

options = optimoptions('fmincon','Display','iter', 'FunctionTolerance', 1e-6, 'ConstraintTolerance', 1e-6, ...
    'Algorithm','active-set', 'OutputFcn', term_fctn); % interior-point, sqp, trust-region-reflective, active-set
options.MaxFunctionEvaluations = 2.5e5;
options.MaxIterations = 1e3;


% Define cost fctn
fun = @(p)cost_fctn(p, A0, b0);

if false
    [p, fval, exitflag, output] = solve_ip(fun, A1, b1, A2, b2, p_init, f0, options);
else
    A2_iter = {};
    for i = 1:nfreq
        A2_iter{end + 1} = zeros(1, size(p_init,1));
    end
    b2_iter = 0.0;

    max_iter = 20;
    tol = 1e3; % Pa
    p = p_init;
    for iter = 1:max_iter
        disp(strcat("Main iteration ", num2str(iter), "..."))

        % Optimize and check if ineq constraints fulfilled
        [p, fval, exitflag, output] = solve_ip(fun, A1_cells, b1, A2_iter, b2_iter, p, f0, options);
        b2_new = abs(timeDomainSum(f0,A2_cells,p));

        unful_ineq = b2_new > (b2 + tol);
        if ~any(unful_ineq)
            break; % successful
        else
            % Add inequalities that are not fulfilled
            for i = 1:nfreq
                A2_iter{i} = [A2_iter{i}; A2_cells{i}(unful_ineq, :)];
            end
            b2_iter = [b2_iter; b2(unful_ineq)];
        end
    end

    if any(unful_ineq)
        disp("Optimization failed");
    else
        disp("Optimization successful");
    end
end

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


function val = cost_fctn(p, A0, b0)

p = getCompVec(p);
% Calc cost fctn value
val = norm(abs(A0 * p) - b0);

end

function [c,ceq] = ineq_const(p, A1_cells, b1, A2_cells, b2, f0)

p = getCompVec(p);

% Transform resulting signals into time domain and sum signals
y1 = timeDomainSum(f0,A1_cells,p);
y2 = timeDomainSum(f0,A2_cells,p);

% Calc equality constraint
ceq = abs(y1) - b1;
% Calc inequality constraint value (if c > 0, then constraint is not satisfied)
c = abs(y2) - b2;

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

function [p, fval, exitflag, output] = solve_ip(fun, A1_cells, b1, A2_cells, b2, p_init, f0, options)

% Update constraints
nonlcon = @(p_init)ineq_const(p_init, A1_cells, b1, A2_cells, b2, f0);

[p_opt, fval, exitflag, output] = fmincon(fun, [real(p_init); imag(p_init)], [], [], [], [], [], [], nonlcon, options);

p = getCompVec(p_opt);

end
