function [A1, A2, b1, b2, A_zero, A_vol] = prepare_opt_vars(A, b, considered_ids, vol_ids, init_ids)

A1 = double(A(considered_ids &  init_ids, :)); % Init points
A2 = double(A(considered_ids & ~init_ids, :)); % Non-init points
b1 = double(b(considered_ids &  init_ids));
b2 = double(b(considered_ids & ~init_ids));

A_zero = double(A(considered_ids & ~vol_ids, :)); % Off-Target Domain
A_vol  = double(A(considered_ids &  vol_ids, :)); % Target Volume

% A1 = A_vol;
% A2 = A_zero;
% b1 = double(b(opt_ids &  obs_ids));
% b2 = double(b(opt_ids & ~obs_ids));

end