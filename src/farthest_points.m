function logical_init_ids = farthest_points(grid_sz, stim_ids, min_dist, hull_ids, cluster_labels)

min_dist_hull = 3;

if length(grid_sz) == 2
    [x, y] = ind2sub(grid_sz, stim_ids);

    % Compute pairwise distances using vectorized operations
    [X1, X2] = meshgrid(x, x);
    [Y1, Y2] = meshgrid(y, y);
    distances = sqrt((X1 - X2).^2 + (Y1 - Y2).^2);
else
    [x, y, z] = ind2sub(grid_sz, stim_ids);

    % Compute pairwise distances using vectorized operations
    [X1, X2] = meshgrid(x, x);
    [Y1, Y2] = meshgrid(y, y);
    [Z1, Z2] = meshgrid(z, z);
    distances = sqrt((X1 - X2).^2 + (Y1 - Y2).^2 + (Z1 - Z2).^2);
end

% Limit the selected points with the surface points
init_ids = [];
selected = false(length(x), 1);

hulls = false(length(x), 1);
hulls(ismember(stim_ids, cat(2, hull_ids{:}))) = true; % Exclude surfaces indices

prev_len = 0;
for i = 1:length(hull_ids) % for each cluster
    selected_i = selected;
    selected_i(cluster_labels ~= i) = true; % Exclude indices from other clusters

    hulls_i = hulls;
    hulls_i(cluster_labels ~= i) = true; % Exclude indices from other clusters

    point_found = false;

    % Compute the minimum distance to the current set of selected points and hulls
    minDistToSelected = min(distances(selected_i | hulls_i, ~selected_i), [], 1);
    minDistToHulls = minDistToSelected;
    
    % Iteratively select the farthest point
    while true        
    
        % Exclude indices where distance to hull below threshold
        minDistToSelected(minDistToHulls < min_dist_hull) = -1.0;

        % Find the index of the farthest point that meets the minimum distance criterion
        [maxDist, nextIndex] = max(minDistToSelected);
        [maxDistHull, ~] = max(minDistToHulls);

        if isempty(maxDist)
            if ~point_found
                % take first index if no point found inside surface
                init_ids(end+1) = find(~selected_i, 1);
            end
            break;
        end
        if maxDist < min_dist(i) || maxDistHull < min_dist_hull
            if ~point_found
                % take the point found inside surface
                nonSelectedIndices = find(~selected_i);
                nextPoint = nonSelectedIndices(nextIndex);
                init_ids(end+1) = nextPoint;
            end
            break;
        end
    
        nonSelectedIndices = find(~selected_i);
        nextPoint = nonSelectedIndices(nextIndex);
    
        % Update the indices and selected array
        init_ids(end+1) = nextPoint;
        selected_i(nextPoint) = true;
        point_found = true;

        % Compute the minimum distance to the current set of selected points and hulls
        minDistToSelected = min(distances(selected_i, ~selected_i), [], 1);
        minDistToHulls = min(distances(hulls_i, ~selected_i), [], 1);
    end

    fprintf(strcat("\nNo of Init ids cluster #", num2str(i), ": ", num2str(length(init_ids) - prev_len)))
    prev_len = length(init_ids);
end
fprintf("\n")

if length(grid_sz) == 2
    init_ids = sub2ind(grid_sz, x(init_ids), y(init_ids));
else
    init_ids = sub2ind(grid_sz, x(init_ids), y(init_ids), z(init_ids));
end

logical_init_ids = false(prod(grid_sz), 1);
logical_init_ids(init_ids) = true;

end
