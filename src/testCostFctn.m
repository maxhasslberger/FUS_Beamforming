function val = testCostFctn(A, b, domain_ids, skull_ids, vol_ids, p_init, init_ids, beta, num_els)

p_init = double(p_init);

% Separate A and b into target region (A1,b1) and other indices (A2,b2)
% Needs to be fixed to account for larger A matrix
% Using cell array for A now

% Prepare opt vars
m = size(A, 2);
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

% Transform resulting signals into time domiain
mv_prod_f1 = A1_cells{1} * p_init(:,1);
ifft_prod_f1 = ifft(mv_prod_f1);

mv_prod_f2 = A1_cells{2} * p_init(:,2);
ifft_prod_f2 = ifft(mv_prod_f2);

% Sum signals in the time domain
ifft_prod = ifft_prod_f2 + ifft_prod_f1;
% Obtain amplitude
y = abs(ifft_prod);

% Compare amplitude of signals with desired pressures

val = y - b1;


% % val = abs(sum(ifft(A1 * p))) - b1;


