% Assuming medium.density is already defined in your workspace

% Step 1: Identify the skull pixels (density > 1000)
skullMask = medium.density > 1000;

% Find the row and column indices of the skull pixels
[rowSkull, colSkull] = find(skullMask);

% Step 2: Use convhull to create a convex hull around the skull pixels
if ~isempty(rowSkull)
    k = convhull(colSkull, rowSkull);
else
    k = [];
end

% Step 3: Create a mask for the convex hull region
convexHullMask = poly2mask(colSkull(k), rowSkull(k), size(medium.density, 1), size(medium.density, 2));

% Step 4: Mark the interior of the convex hull
interiorMask = convexHullMask & ~skullMask;

% Find the linear indices of the interior and skull pixels
interiorIndices = find(interiorMask);
skullIndices = find(skullMask);

% Convert linear indices to subscripts (row and column indices)
[rowInterior, colInterior] = ind2sub(size(medium.density), interiorIndices);
[rowSkull, colSkull] = ind2sub(size(medium.density), skullIndices);

% Plot the density matrix
figure;
imagesc(medium.density);
colormap(gray);
colorbar;
hold on;

% Draw the convex hull around the skull
plot(colSkull(k), rowSkull(k), 'g-', 'LineWidth', 2);

% Mark the skull region on the plot
plot(colSkull, rowSkull, 'b.', 'MarkerSize', 15);

% Mark the interior region on the plot
plot(colInterior, rowInterior, 'r.', 'MarkerSize', 15);

% Add labels and title for better understanding
xlabel('Column Index');
ylabel('Row Index');
title('Density Matrix with Skull and Interior Region Encapsulated by Convex Hull');

% Display the results
disp('Row indices of pixels inside the convex hull:');
disp(rowInterior);
disp('Column indices of pixels inside the convex hull:');
disp(colInterior);

disp('Row indices of skull pixels:');
disp(rowSkull);
disp('Column indices of skull pixels:');
disp(colSkull);
