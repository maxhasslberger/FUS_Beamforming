function init_ids = get_init_ids(kgrid, lambda, b_mask)

stim_ids = find(b_mask);
min_dist = 2 * lambda / kgrid.dx; % m -> pix; min distance between init points

% Find n farthest points
init_ids = farthest_points(size(b_mask), stim_ids, min_dist);

% Plot selected points
b_mask(init_ids) = 5;
plot_results(kgrid, [], b_mask, 'Initial Points', [], [], [1 1 1], size(b_mask), 1, false, [])

end
