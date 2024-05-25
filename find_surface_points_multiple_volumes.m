function [surface_indices, labels] = find_surface_points_multiple_volumes(epsilon, x, y, z)
    % Input:
    % x - Row vector of x coordinates of the points
    % y - Column vector of y coordinates of the points
    % z - (Optional) Column vector of z coordinates of the points
    % epsilon - The maximum distance between two points for one to be considered in the same cluster

    % Check if z is provided (3D case)
    if nargin < 4 || isempty(z)
        is3D = false;
        z = [];
    else
        is3D = true;
    end

    % Ensure x, y, and z are column vectors
    x = x(:);
    y = y(:);
    if is3D
        z = z(:);
    end

    % Combine x, y, and z (if 3D) into a single matrix
    if is3D
        data = [x, y, z];
    else
        data = [x, y];
    end
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
            if is3D %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Clean up and test 3D case (sliceViewer)
                % Compute the convex hull for the 3D cluster
                hull_indices = convhull(cluster_points(:, 1), cluster_points(:, 2), cluster_points(:, 3));
                hull_indices = unique(hull_indices); % only memorize ids contained in convex hull
%                 k = convhulln([cluster_points(:, 1), cluster_points(:, 2), cluster_points(:, 3)]);
            else
                % Compute the convex hull for the 2D cluster
                hull_indices = convhull(cluster_points(:, 1), cluster_points(:, 2));
                hull_indices = hull_indices(1:end-1); % Remove the duplicate closing point
            end
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
