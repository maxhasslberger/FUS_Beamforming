function init_ids = get_init_ids(kgrid, lambda, b_mask)

stim_ids = find(b_mask);
min_dist = 1 * lambda / kgrid.dx; % m -> pix; min distance between init points

% Obtain hull points of sonicated volumes
epsilon = 2;
[x, y] = ind2sub(size(b_mask), stim_ids);
[hull_ids, cluster_labels] = get_hulls_multiple_volumes(epsilon, x, y);

% Find n farthest points
init_ids = farthest_points(size(b_mask), stim_ids, min_dist, hull_ids, cluster_labels);

% Plot selected points
b_mask(init_ids) = 5;
for i = 1:length(hull_ids)
    [x, y] = ind2sub(size(b_mask), stim_ids);
    plot_hull_ids = sub2ind(size(b_mask), x(hull_ids{i}), y(hull_ids{i}));
    b_mask(plot_hull_ids) = 3;
end

plot_results(kgrid, [], b_mask, 'Initial Points', [], [], [1 1 1], size(b_mask), 1, false, [])

end
