function [A1, A2] = prepare_opt_vars(A, considered_ids, vol_ids, init_ids)

A1 = double(A(considered_ids &  init_ids, :)); % Init points
A2 = double(A(considered_ids & ~init_ids, :)); % Non-init points