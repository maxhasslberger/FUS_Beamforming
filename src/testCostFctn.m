function [val,num_failed_constraints] = testCostFctn(A, b, domain_ids, skull_ids, vol_ids, p_init, init_ids, beta, num_els)
% Test cost function and inequality constraint logic for multiple frequencies using initial solution
% Does not include fmincon optimization so that we are only working with one iteration for testing

p_init = double(p_init);

% Separate A and b into target region (A1,b1) and other indices (A2,b2)
% Needs to be fixed to account for larger A matrix
% Using cell array for A now


m = size(A, 2);

% Prepare opt vars
A1_cells = {};
A2_cells = {};
for i = 1:m
    A_temp = A{i};
    [A1, A2] = prepare_opt_A(A_temp, domain_ids | skull_ids, vol_ids, init_ids);
    % A1 = double(A_temp(domain_ids | skull_ids &  init_ids, :)); % Init points
    % A2 = double(A_temp(domain_ids | skull_ids & ~init_ids, :)); % Non-init points
    A1_cells{end + 1} = A1;
    A2_cells{end + 1} = A2;
end
[b1, b2] = prepare_opt_b(b, domain_ids | skull_ids, vol_ids, init_ids);


% Add regularization
for i = 1:m
    [A1_cells{i}, b1] = add_L2_reg(A1_cells{i}, b1, beta(1));
    [A2_cells{i}, b2] = add_ineq(A2_cells{i}, b2, length(p_init));
end

% Transform resulting signals into time domain and sum signals
y1 = timeDomainSum(m,A1_cells,p_init);
y2 = timeDomainSum(m,A2_cells,p_init);

% Calc cost fctn value
val = norm(abs(y1) - b1);

% Calc inequality constraint value (if c > 0, then constraint is not satisfied)
c = abs(y2) - b2;
num_failed_constraints = length(find(c > 0));



% Old code:
% for i = 1:m
%     y_temp = ifft(A1_cells{i} * p_init(:,i));
%     y = [y + y_temp];
% end

% Compare amplitude of signals with desired pressures


end


% function y = timeDomainSum(m,A_cells,p)
%     n = size(A_cells{1},1)
%     y = zeros(n,1);
%     for i = 1:m
%         y_temp = ifft(A_cells{i} * p(:,i));
%         y = [y + y_temp];
%     end
% end




