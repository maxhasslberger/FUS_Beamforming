clear;

f0 = 500e3; % Hz - transducer frequency
n_dim = 3;
dx_factor = 1;
[kgrid, medium, ppp] = init_grid_medium(f0, 'n_dim', n_dim, 'dx_factor', 1 / dx_factor);

% Init
karray_t = kWaveArray();
element_pos = load(fullfile('Array_Positions\std.mat')).ElementPosition';

% Add one array element after another
for i = 1:length(element_pos)
%     karray_t.addCustomElement(round(element_pos(:, i), 3), 0, 2, char("el_" + string(i)));
    karray_t.addDiscElement(round(element_pos(:, i), 3), 1e-10, ones(1, 3));
%     mask = karray_t.getArrayBinaryMask(kgrid);
%     voxelPlot(double(mask))
end

% TODO: Rotation and translation considered before added to karray class

% % Translation and Rotation of every single element...
% translation = [0, 0, 0] * 1e-3;
% rotation = [0, 90, 0]; % only integers of 90° supported
% karray_t.setArrayPosition(translation, rotation);

mask = karray_t.getArrayBinaryMask(kgrid);
voxelPlot(double(mask))
