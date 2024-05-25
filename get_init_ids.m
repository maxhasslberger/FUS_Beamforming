function init_ids = get_init_ids(kgrid, lambda, b_mask)

stim_ids = find(b_mask);
min_dist = 2 * lambda / kgrid.dx; % m -> pix; min distance between init points

% Obtain surface points of sonicated volumes
epsilon = 2;
[x, y] = ind2sub(size(b_mask), stim_ids);
[surface_ids, cluster_labels] = find_surface_points_multiple_volumes(epsilon, x, y);

% Find n farthest points
init_ids = farthest_points(size(b_mask), stim_ids, min_dist, surface_ids, cluster_labels);

% Plot selected points
b_mask(init_ids) = 5;
for i = 1:length(surface_ids)
    [x, y] = ind2sub(size(b_mask), stim_ids);
    b_mask(x(surface_ids{i}), y(surface_ids{i})) = 3; %%%%%%%%%%%%%%%%%%%%% Visualize surfaces
end

plot_results(kgrid, [], b_mask, 'Initial Points', [], [], [1 1 1], size(b_mask), 1, false, [])

end
