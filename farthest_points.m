function indices = farthest_points(x, y, n_points)

x = x(:);
y = y(:);

% Compute pairwise distances using vectorized operations
[X1, X2] = meshgrid(x, x);
[Y1, Y2] = meshgrid(y, y);
distances = sqrt((X1 - X2).^2 + (Y1 - Y2).^2);

% Initialize the selected points with the first point
indices = zeros(1, n_points);
indices(1) = 1;
selected = false(length(x), 1);
selected(1) = true;

% Iteratively select the farthest point
for k = 2:n_points
    % Compute the minimum distance to the current set of selected points
    minDistToSelected = min(distances(selected, ~selected), [], 1);

    % Find the index of the farthest point
    [~, nextIndex] = max(minDistToSelected);
    nonSelectedIndices = find(~selected);
    nextPoint = nonSelectedIndices(nextIndex);

    % Update the indices and selected array
    indices(k) = nextPoint;
    selected(nextPoint) = true;
end

end
