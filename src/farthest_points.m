function logical_init_ids = farthest_points(grid_sz, stim_ids, min_dist, hull_ids, cluster_labels)

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
selected(ismember(stim_ids, cat(2, hull_ids{:}))) = true; % Exclude surfaces indices

for i = 1:length(hull_ids)
    selected_i = selected;
    selected_i(cluster_labels ~= i) = true; % Exclude indices from other clusters

    point_found = false;
    
    % Iteratively select the farthest point
    while true
        % Compute the minimum distance to the current set of selected points
        minDistToSelected = min(distances(selected_i, ~selected_i), [], 1);
    
        % Find the index of the farthest point that meets the minimum distance criterion
        [maxDist, nextIndex] = max(minDistToSelected);
        if isempty(maxDist)
            if ~point_found
                % take first index if no point found inside surface
                init_ids(end+1) = find(~selected_i, 1);
            end
            break;
        end
        if maxDist < min_dist
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
    end
end

if length(grid_sz) == 2
    init_ids = sub2ind(grid_sz, x(init_ids), y(init_ids));
else
    init_ids = sub2ind(grid_sz, x(init_ids), y(init_ids), z(init_ids));
end

logical_init_ids = false(prod(grid_sz), 1);
logical_init_ids(init_ids) = true;

end
