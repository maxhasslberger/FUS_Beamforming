function init_ids = get_farthest_ids(kgrid, lambda, b_mask)

stim_ids = find(b_mask);
min_dist = 2 * lambda / kgrid.dx; % m -> pix; min distance between init points

% Find n farthest points
[stim_ids_r, stim_ids_c] = ind2sub(size(b_mask), stim_ids);
rc_ids = farthest_points(stim_ids_c, stim_ids_r, min_dist);
init_ids = sub2ind(size(b_mask), stim_ids_r(rc_ids), stim_ids_c(rc_ids));

% Plot selected points
b_mask(init_ids) = 5;
plot_results(kgrid, [], b_mask, 'Initial Points', [], [], [1 1 1], size(b_mask), 1, false, [])

end
