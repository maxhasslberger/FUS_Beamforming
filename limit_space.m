function logical_opt_ids = limit_space(sound_speed)

dim = numel(size(sound_speed));
c0 = min(sound_speed(:));

if max(sound_speed(:)) == c0 % homogeneous medium
    logical_opt_ids = true(numel(sound_speed), 1);
    return;
end

skullMask = sound_speed > c0;

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

    % Mark the skull region on the plot
    plot(colSkull, rowSkull, 'b.', 'MarkerSize', 15);
    
    % Add labels and title for better understanding
    xlabel('Column Index');
    ylabel('Row Index');
    title('Density Matrix with Skull and Interior Region Encapsulated by Convex Hull');
    
%     % Display the results
%     disp('Row indices of pixels inside the convex hull:');
%     disp(row_opt);
%     disp('Column indices of pixels inside the convex hull:');
%     disp(col_opt);
%     
%     disp('Row indices of skull pixels:');
%     disp(rowSkull);
%     disp('Column indices of skull pixels:');
%     disp(colSkull);
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
        k = [];
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
    logical_opt_ids = false(sum(size(sound_speed)), 1);
    logical_opt_ids(opt_ids) = true;
    
    % Visualize
    % Find the linear indices of the interior and skull pixels
    skullIndices = find(skullMask);
    
    % Convert linear indices to subscripts (row, column, and depth indices)
    [row_opt, col_opt, depth_opt] = ind2sub(size(sound_speed), opt_ids);
    [rowSkull, colSkull, depthSkull] = ind2sub(size(sound_speed), skullIndices);
    sliceToShow = 5; % Choose the slice to display
    
    % Plot the density matrix (sliceToShow)
    %figure;
    %imagesc(sound_speed(:,:,sliceToShow));
    %colormap(gray);
    %hold on;
    
    % Draw the convex hull around the skull on the slice
    if false%~isempty(k)
        plot3(colSkull(k), rowSkull(k), 'g-', 'LineWidth', 2);
    end

    figure;
    sliceViewer(convexHullMask + skullMask, 'SliceDirection','X')
    
    % Mark the interior region on the plot
    %plot3(col_opt, row_opt, 'r.', 'MarkerSize', 15);

    % Mark the skull region on the plot
    %plot3(colSkull, rowSkull, 'b.', 'MarkerSize', 15);
    
    % Add labels and title for better understanding
    %xlabel('Column Index');
    %ylabel('Row Index');
    %title('Density Matrix (Slice) with Skull and Interior Regions Encapsulated by Convex Hull');
    
%     % Display the results
%     disp('Row indices of pixels inside the convex hull:');
%     disp(row_opt);
%     disp('Column indices of pixels inside the convex hull:');
%     disp(col_opt);
%     disp('Depth indices of pixels inside the convex hull:');
%     disp(depth_opt);
%     
%     disp('Row indices of skull pixels:');
%     disp(rowSkull);
%     disp('Column indices of skull pixels:');
%     disp(colSkull);
%     disp('Depth indices of skull pixels:');
%     disp(depthSkull);
end
    

end

% Helper function to check if points are inside a convex hull
function inside = isInsideConvexHull(dt, points)
    [~, in] = pointLocation(dt, points);
    inside = ~isnan(in);
end

function inside = pointInTetrahedron(points, tetra)
    % Create matrices for the point-in-tetrahedron test
    A = tetra(1, :);
    B = tetra(2, :);
    C = tetra(3, :);
    D = tetra(4, :);
    
    % Compute vectors
    AB = B - A;
    AC = C - A;
    AD = D - A;
    AP = points - A;
    
    % Compute dot products
    AB_AB = dot(AB, AB, 2);
    AB_AC = dot(AB, AC, 2);
    AB_AD = dot(AB, AD, 2);
    AB_AP = dot(AB, AP, 2);
    AC_AC = dot(AC, AC, 2);
    AC_AD = dot(AC, AD, 2);
    AC_AP = dot(AC, AP, 2);
    AD_AD = dot(AD, AD, 2);
    AD_AP = dot(AD, AP, 2);
    
    % Compute volumes
    denom = AB_AB * (AC_AC * AD_AD - AC_AD * AC_AD) + AB_AC * (AB_AC * AD_AD - AB_AD * AC_AD) + AB_AD * (AB_AC * AC_AD - AB_AC * AD_AD);
    v = (AC_AC * AD_AD - AC_AD * AC_AD) * AB_AP + (AB_AC * AD_AD - AB_AD * AC_AD) * AC_AP + (AB_AC * AC_AD - AB_AD * AC_AC) * AD_AP;
    w = (AB_AB * AD_AD - AB_AD * AB_AD) * AC_AP + (AB_AB * AC_AD - AB_AC * AB_AD) * AD_AP + (AB_AC * AB_AD - AB_AB * AC_AD) * AB_AP;
    u = (AB_AB * (AC_AC * AD_AP - AC_AD * AC_AP) + AB_AC * (AB_AC * AD_AP - AB_AD * AC_AP) + AB_AD * (AB_AC * AC_AP - AB_AB * AD_AP));
    
    % Check if point is inside the tetrahedron
    inside = u >= 0 & v >= 0 & w >= 0 & (u + v + w) <= denom;
end
