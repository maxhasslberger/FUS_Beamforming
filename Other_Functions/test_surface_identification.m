% Example usage for 2D case
x2 = [1, 2, 3, 4, 5, 3, 3, 10, 11, 12, 12, 11];
y2 = [5, 4, 3, 2, 1, 4, 2, 10, 11, 10, 9, 8];
epsilon2 = 2; % Adjust epsilon based on the density of your data

[surface_indices2, labels2] = find_surface_points_multiple_volumes(epsilon2, x2, y2);

disp('2D Surface indices for each volume:');
for i = 1:length(surface_indices2)
    disp(['Volume ' num2str(i) ':']);
    disp(surface_indices2{i});
end

% Example usage for 3D case
x3 = [1, 2, 3, 4, 5, 3, 3, 10, 11, 12, 12, 11];
y3 = [5, 4, 3, 2, 1, 4, 2, 10, 11, 10, 9, 8];
z3 = [1, 1, 1, 1, 1, 2, 2, 3, 3, 3, 4, 4];
epsilon3 = 2; % Adjust epsilon based on the density of your data

[surface_indices3, labels3] = find_surface_points_multiple_volumes(epsilon3, x3, y3, z3);

disp('3D Surface indices for each volume:');
for i = 1:length(surface_indices3)
    disp(['Volume ' num2str(i) ':']);
    disp(surface_indices3{i});
end
