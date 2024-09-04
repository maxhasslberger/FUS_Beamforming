function [init_ids, hull_ids, b_mask_plot] = get_init_ids(kgrid, min_dist, b_mask, force_pressures)

n_dim = length(size(kgrid.k));
min_dist = min_dist(force_pressures) / kgrid.dx; % m -> pix

% Label clusters and identify forced stim_ids
stim_ids = [];
cluster_labels = [];

for i = 1:length(force_pressures)
    if n_dim == 2
        new_ids = find(b_mask(:, :, force_pressures(i)));
    else
        new_ids = find(b_mask(:, :, :, force_pressures(i)));
    end

    stim_ids = [stim_ids; new_ids];
    cluster_labels = [cluster_labels; i * ones(length(new_ids), 1)];
end

% Obtain hull points of sonicated volumes
hull_ids = get_hulls(cluster_labels, size(kgrid.k), stim_ids);

% Find n farthest points
init_ids = farthest_points(size(kgrid.k), stim_ids, min_dist, hull_ids, cluster_labels);

% Plot selected points
b_mask_plot = zeros(size(kgrid.k));
b_mask_plot(init_ids) = 1;
for i = 1:length(hull_ids)
%     [x, y] = ind2sub(size(b_mask_plot), stim_ids);
%     plot_hull_ids = sub2ind(size(b_mask_plot), x(hull_ids{i}), y(hull_ids{i}));
%     b_mask_plot(plot_hull_ids) = 0.5;
    b_mask_plot(hull_ids{i}) = 0.5;
end

% figure;
% sliceViewer(b_mask_plot, 'SliceDirection','X')

end
