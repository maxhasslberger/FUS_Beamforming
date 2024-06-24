function [logical_opt_ids, logical_skull_ids] = limit_space(sound_speed)

dim = numel(size(sound_speed));
c0 = min(sound_speed(:));

if max(sound_speed(:)) == c0 % homogeneous medium
    logical_opt_ids = true(numel(sound_speed), 1);
    return;
end

skullMask = sound_speed > c0;
logical_skull_ids = reshape(skullMask, [], 1);
skull_ids = find(skullMask);

if dim == 2
    % Find the row and column indices of the skull pixels
    [rowSkull, colSkull] = find(skullMask);

    % Create a convex hull around the skull pixels
    if ~isempty(rowSkull)
        k = convhull(colSkull, rowSkull);
    else
        k = [];
    end
    
    % Create a mask and obtain indices for the convex hull region
    convexHullMask = poly2mask(colSkull(k), rowSkull(k), size(sound_speed, 1), size(sound_speed, 2));
    opt_ids = find(convexHullMask);
    opt_ids = setdiff(opt_ids, skull_ids);
    logical_opt_ids = false(numel(sound_speed), 1);
    logical_opt_ids(opt_ids) = true;
    
    % Visualize
    [row_opt, col_opt] = ind2sub(size(sound_speed), opt_ids);
    
    figure;
    imagesc(sound_speed);
    colormap(gray);
    hold on;
    
    % Draw the convex hull around the skull
    if ~isempty(k)
        plot(colSkull(k), rowSkull(k), 'g-', 'LineWidth', 2);
    end
    
    % Mark the interior region on the plot
    plot(col_opt, row_opt, 'r.', 'MarkerSize', 15);
    
    % Add labels and title for better understanding
    xlabel('Column Index');
    ylabel('Row Index');
    title('Density Matrix with Skull and Interior Region Encapsulated by Convex Hull');
else
    % Find the row, column, and depth indices of the skull pixels
    [rowSkull, colSkull, depthSkull] = ind2sub(size(sound_speed), find(skullMask));
    
    % Create a convex hull around the skull pixels
    if ~isempty(rowSkull)
        % Create a point cloud of skull pixels
        skullPointCloud = [rowSkull, colSkull, depthSkull];
        
        % Compute the convex hull in 3D
        k = convhull(rowSkull, colSkull, depthSkull);
        k = unique(k(:)); % only memorize ids contained in convex hull
        hullPoints = skullPointCloud(k, :); % Points forming the convex hull
    else
        hullPoints = [];
    end
    
    % Create a mask and obtain indices for the convex hull region
    convexHullMask = false(size(sound_speed));
    if ~isempty(hullPoints)
        % Create a triangulation object
        dt = delaunayTriangulation(hullPoints);
        
        % Get all points in the 3D space of the matrix
        [X, Y, Z] = ndgrid(1:size(sound_speed, 1), 1:size(sound_speed, 2), 1:size(sound_speed, 3));
        testPoints = [X(:), Y(:), Z(:)];
        
        % Check which points are inside the convex hull
        insideHull = isInsideConvexHull(dt, testPoints);
        
        % Reshape to the original size of the 3D matrix
        convexHullMask = reshape(insideHull(:, 1), size(sound_speed));
    end
    opt_ids = find(convexHullMask);
    opt_ids = setdiff(opt_ids, skull_ids);
    logical_opt_ids = false(numel(sound_speed), 1);
    logical_opt_ids(opt_ids) = true;
    
    % Visualize the masked region
    convexHullMask = zeros(size(sound_speed));
    convexHullMask(opt_ids) = 1;
    figure;
    sliceViewer(convexHullMask, 'SliceDirection','X')
    
end
    

end

% Helper function to check if points are inside a convex hull
function inside = isInsideConvexHull(dt, points)
    [~, in] = pointLocation(dt, points);
    inside = ~isnan(in);
end

