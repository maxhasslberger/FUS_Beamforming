function [hull_ids, labels] = get_hulls_multiple_volumes(epsilon, grid_sz, stim_ids)

% Obtain coordinates of stim_ids
if length(grid_sz) == 2 % dim == 2
    [x, y] = ind2sub(grid_sz, stim_ids);
    data = [x, y];
else
    [x, y, z] = ind2sub(grid_sz, stim_ids);
    data = [x, y, z];
end

original_indices = (1:length(x))';
numPoints = length(x);

% Initialize labels
labels = zeros(numPoints, 1);
cluster_id = 0;

% Custom distance-based clustering
for i = 1:numPoints
    if labels(i) == 0
        cluster_id = cluster_id + 1;
        labels = assign_cluster(data, labels, i, cluster_id, epsilon);
    end
end

% Get unique cluster labels
unique_labels = unique(labels);

% Initialize a cell array to hold surface indices of each cluster
hull_ids = cell(length(unique_labels), 1);

% Loop over each cluster and compute the convex hull
for i = 1:length(unique_labels)
    cluster_label = unique_labels(i);
    cluster_indices = original_indices(labels == cluster_label);
    cluster_points = data(labels == cluster_label, :);

    if size(cluster_points, 1) >= 3
        k = convhull(cluster_points);
        k = unique(k(:)); % only memorize unique ids contained in convex hull
        hull_ids{i} = cluster_indices(k)';
    else
        % If less than 3 points, all points are part of the convex hull
        hull_ids{i} = cluster_indices';
    end

    % Sort the surface indices according to the original indices
    hull_ids{i} = sort(hull_ids{i});
end

end

function labels = assign_cluster(data, labels, point_index, cluster_id, epsilon)

% Custom function to assign cluster id based on epsilon distance
labels(point_index) = cluster_id;
for j = 1:length(data)
    if labels(j) == 0
        distance = norm(data(point_index, :) - data(j, :));
        if distance <= epsilon
            labels = assign_cluster(data, labels, j, cluster_id, epsilon);
        end
    end
end

end
