function [karray_t, el2mask_ids] = create_transducer(kgrid, t_name, t_pos, t_rot)

% Init
element_pos = load(fullfile("Array_Positions", t_name + ".mat")).ElementPosition'; % flat transducer array centered at [0, 0, 0] along the xy-plane
karray_t = kWaveArray();

n_tr_elements = size(element_pos, 2);
elementAll_pos = nan(3, n_tr_elements * size(t_pos, 2));
for tx = 1:size(t_pos, 2) % for each transducer

    tx_pos = t_pos(:, tx);
    tx_rot = t_rot(:, tx);

    % Shift and align transducers 
    tx_transf = getAffineMatrix(tx_pos, tx_rot);
    
    elementx_pos = tx_transf(1:3, 1:3) * element_pos;
    elementAll_pos(:, (tx - 1) * n_tr_elements + 1 : tx * n_tr_elements) = elementx_pos + repmat(tx_pos, 1, size(element_pos, 2));
end

[~, el2mask_ids] = sortrows(elementAll_pos'); % refer element indices to mask
    
% Add one array element after another
for i = 1:length(elementAll_pos)
%     karray_t.addCustomElement(round(element_pos(:, i), 3), 0, 2, char("el_" + string(i)));
    karray_t.addDiscElement(round(elementAll_pos(:, i), 3), 1e-10, zeros(1, 3));
%     mask = karray_t.getArrayBinaryMask(kgrid);
%     voxelPlot(double(mask))
end
% mask = karray_t.getArrayBinaryMask(kgrid);
% voxelPlot(double(mask))

end
