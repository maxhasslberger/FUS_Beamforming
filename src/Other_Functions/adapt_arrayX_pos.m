% Function to adapt the spacing for transducer arrayX to read the mask in karray class

element_pos_orig = load(fullfile('..\..\Array_Positions\std_orig.mat')).ElementPosition';
element_pos = round(element_pos_orig, 3);
shiftx = zeros(1, length(element_pos));
shifty = zeros(1, length(element_pos));

shiftx(element_pos(1, :) > 0) = -1e-3;
shifty(element_pos(2, :) > 0) = -1e-3;

element_pos(1, :) = element_pos(1, :) + shiftx;
element_pos(2, :) = element_pos(2, :) + shifty;

scatter(element_pos_orig(1,:), element_pos_orig(2,:), 'filled')
figure
scatter(element_pos(1,:), element_pos(2,:), 'filled')

ElementPosition = element_pos';
save(fullfile('..\..\Array_Positions\std.mat'), "ElementPosition")
