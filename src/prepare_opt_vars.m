function [A1, A2, b1, b2] = prepare_opt_vars(A, b, considered_ids, vol_ids, init_ids, ineq_active)

% Equality constraints
eq_ids = considered_ids &  init_ids;
A1 = double(A(eq_ids, :));
b1 = double(b(eq_ids));

% Inequality constraints
if ineq_active
    ineq_ids = considered_ids & ~init_ids;
else
    ineq_ids = vol_ids & ~init_ids; % Still consider limited volumes
end

A2 = double(A(ineq_ids, :));
b2 = double(b(ineq_ids));

% A_zero = double(A(considered_ids & ~vol_ids, :)); % Off-Target Domain
% A_vol  = double(A(considered_ids &  vol_ids, :)); % Target Volume

% A1 = A_vol;
% A2 = A_zero;
% b1 = double(b(opt_ids &  obs_ids));
% b2 = double(b(opt_ids & ~obs_ids));

end