function indices = farthest_points(x, y, min_dist) % TODO: 3D

x = x(:);
y = y(:);

% Compute pairwise distances using vectorized operations
[X1, X2] = meshgrid(x, x);
[Y1, Y2] = meshgrid(y, y);
distances = sqrt((X1 - X2).^2 + (Y1 - Y2).^2);

% Initialize the selected points with the first point
indices = [];
indices(1) = 1;
selected = false(length(x), 1);
selected(1) = true;

% Iteratively select the farthest point
while true
    % Compute the minimum distance to the current set of selected points
    minDistToSelected = min(distances(selected, ~selected), [], 1);

    % Find the index of the farthest point that meets the minimum distance criterion
    [maxDist, nextIndex] = max(minDistToSelected);
    if maxDist < min_dist
        break;
    end

    nonSelectedIndices = find(~selected);
    nextPoint = nonSelectedIndices(nextIndex);

    % Update the indices and selected array
    indices(end+1) = nextPoint;
    selected(nextPoint) = true;
end

end
