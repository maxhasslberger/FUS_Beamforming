function init_ids = farthest_points(grid_sz, stim_ids, min_dist) % TODO: 3D

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

% Initialize the selected points with the first point
init_ids = [];
init_ids(1) = 1;
selected = false(length(x), 1);
selected(1) = true;

% Iteratively select the farthest point
while true
    % Compute the minimum distance to the current set of selected points
    minDistToSelected = min(distances(selected, ~selected), [], 1);

    % Find the index of the farthest point that meets the minimum distance criterion
    [maxDist, nextIndex] = max(minDistToSelected);
    if isempty(maxDist)
        break;
    end
    if maxDist < min_dist
        break;
    end

    nonSelectedIndices = find(~selected);
    nextPoint = nonSelectedIndices(nextIndex);

    % Update the indices and selected array
    init_ids(end+1) = nextPoint;
    selected(nextPoint) = true;
end

if length(grid_sz) == 2
    init_ids = sub2ind(grid_sz, x(init_ids), y(init_ids));
else
    init_ids = sub2ind(grid_sz, x(init_ids), y(init_ids), z(init_ids));
end

end
