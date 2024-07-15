function [b1, b2] = prepare_opt_vars(b, considered_ids, vol_ids, init_ids)

b1 = double(b(considered_ids &  init_ids));
b2 = double(b(considered_ids & ~init_ids));
