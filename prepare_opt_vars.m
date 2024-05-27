function [A1, A2, b1, b2, b_zero, b_vol, p_init] = prepare_opt_vars(A, b, p_init, opt_ids, obs_ids, init_ids)

A1 = double(A(opt_ids &  init_ids, :));
A2 = double(A(opt_ids & ~init_ids, :));

b1 = double(b(opt_ids &  init_ids));
b2 = double(b(opt_ids & ~init_ids));
b_zero = double(b(opt_ids & ~obs_ids));
b_vol  = double(b(opt_ids &  obs_ids));

p_init = double(p_init);

end