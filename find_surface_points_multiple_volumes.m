function surface_indices = find_surface_points_multiple_volumes(x, y, epsilon)
    % Input:
    % x - Row vector of x coordinates of the points
    % y - Column vector of y coordinates of the points
    % epsilon - The maximum distance between two points for one to be considered in the same cluster

    % Ensure x and y are column vectors
    x = x(:);
    y = y(:);

    % Combine x and y into a single matrix
    data = [x, y];
    original_indices = (1:length(x))';

    % Number of points
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
    surface_indices = cell(length(unique_labels), 1);

    % Loop over each cluster and compute the convex hull
    for i = 1:length(unique_labels)
        cluster_label = unique_labels(i);
        cluster_indices = original_indices(labels == cluster_label);
        cluster_points = data(labels == cluster_label, :);

        if size(cluster_points, 1) >= 3
            % Compute the convex hull for the cluster
            hull_indices = convhull(cluster_points(:, 1), cluster_points(:, 2));
            hull_indices = hull_indices(1:end-1); % Remove the duplicate closing point
            surface_indices{i} = cluster_indices(hull_indices)';
        else
            % If less than 3 points, all points are part of the convex hull
            surface_indices{i} = cluster_indices';
        end

        % Sort the surface indices according to the original indices
        surface_indices{i} = sort(surface_indices{i});
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
