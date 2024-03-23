clear;

f0 = 500e3; % Hz - transducer frequency
n_dim = 3;
dx_factor = 1;
[kgrid, medium, ppp] = init_grid_medium(f0, 'n_dim', n_dim, 'dx_factor', 1 / dx_factor);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Init
element_pos = load(fullfile('Array_Positions\std.mat')).ElementPosition'; % flat transducer array centered along [0, 0, 0] on the xy-plane
t_pos = [-45, 0; 0, -45; 0, 0] * 1e-3; % m
t_rot = [0, 90; 90, 0; 0, 0]; % deg

karray_t = kWaveArray();

for tx = 1:size(t_pos, 2) % for each transducer

    tx_pos = t_pos(:, tx);
    tx_rot = t_rot(:, tx);

    % Shift and align transducers 
    tx_transf = getAffineMatrix(tx_pos, tx_rot);
    
    elementx_pos = tx_transf(1:3, 1:3) * element_pos;
    elementx_pos = elementx_pos + repmat(tx_pos, 1, size(element_pos, 2));
    
    % Add one array element after another
    for i = 1:length(elementx_pos)
%         karray_t.addCustomElement(round(element_pos(:, i), 3), 0, 2, char("el_" + string(i)));
        karray_t.addDiscElement(round(elementx_pos(:, i), 3), 1e-10, ones(1, 3)); % Rounding error can introduce error in spacing!
%         mask = karray_t.getArrayBinaryMask(kgrid);
%         voxelPlot(double(mask))
    end
end

% TODO: Rotation and translation considered before added to karray class

mask = karray_t.getArrayBinaryMask(kgrid);
voxelPlot(double(mask))
