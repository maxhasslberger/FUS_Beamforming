function [A1, A2, b1, b2, A_zero, A_vol, gamma] = prepare_opt_vars(A, b, opt_ids, obs_ids, init_ids)

A1 = double(A(opt_ids &  init_ids, :));
A2 = double(A(opt_ids & ~init_ids, :));
b1 = double(b(opt_ids &  init_ids));
b2 = double(b(opt_ids & ~init_ids));

A_zero = double(A(opt_ids & ~obs_ids, :));
A_vol  = double(A(opt_ids &  obs_ids, :));

gamma = ones(length(b1), 1);

% A1 = A_vol;
% A2 = A_zero;
% b1 = double(b(opt_ids &  obs_ids));
% b2 = double(b(opt_ids & ~obs_ids));

end