function init_ids = get_farthest_ids(b_mask, n_init_points)

stim_ids = find(b_mask);

% Find n farthest points
[stim_ids_r, stim_ids_c] = ind2sub(size(b_mask), stim_ids);
rc_ids = farthest_points(stim_ids_c, stim_ids_r, n_init_points);
init_ids = sub2ind(size(b_mask), stim_ids_r(rc_ids), stim_ids_c(rc_ids));

b_mask(init_ids) = 5;
kgrid.dim = 2;
kgrid.dx = 1e-3;
kgrid.dy = 1e-3;
plot_results(kgrid, [], b_mask, 'Initial Points', [], [], [1 1 1], size(b_mask), 1, false, [])

end
