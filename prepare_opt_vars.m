function [A1, A2, b1, b2, p_init] = prepare_opt_vars(A, b, p_init, opt_ids, init_ids)

A1 = double(A(opt_ids &  init_ids, :));
A2 = double(A(opt_ids & ~init_ids, :));

b1 = double(b(opt_ids &  init_ids));
b2 = double(b(opt_ids & ~init_ids));

p_init = double(p_init);

end