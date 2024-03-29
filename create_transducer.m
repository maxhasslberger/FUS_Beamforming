function [karray_t, el2mask_ids] = create_transducer(kgrid, t_name, t_pos, t_rot)

% Init
dim = 2;
element_pos = load(fullfile("Array_Positions", t_name + ".mat")).ElementPosition'; % flat transducer array centered at [0, 0, 0] along the xy-plane
% element_pos = element_pos(1:dim, [1, 13, 26, 39, 42, 55, 67, 140, 210])%2D
% element_pos = element_pos(1:dim, [1, 67, 140, 210])%2D
element_pos = element_pos(1:dim, 1:140);%2D
karray_t = kWaveArray();

n_tr_elements = size(element_pos, 2);
elementAll_pos_orig = nan(dim, n_tr_elements * size(t_pos, 2));
for tx = 1:size(t_pos, 2) % for each transducer

    tx_pos = t_pos(:, tx);
    tx_rot = t_rot(:, tx);

    % Shift and align transducers 
    tx_transf = getAffineMatrix(tx_pos, tx_rot);
    
    elementx_pos = tx_transf(1:dim, 1:dim) * element_pos;
    elementAll_pos_orig(:, (tx - 1) * n_tr_elements + 1 : tx * n_tr_elements) = elementx_pos + repmat(tx_pos, 1, size(element_pos, 2));
end

[elementAll_pos, el2mask_ids] = sortrows(elementAll_pos_orig', 2); % refer element indices to mask -> getDistributedSourceSignal
elementAll_pos = elementAll_pos'; % TODO: Fix indexing!
% elementAll_pos = elementAll_pos_orig;
elementAll_pos

% TODO: Real mask with all elements (karray) + mask with center elements only
    
% Add one array element after another
for i = 1:length(elementAll_pos)
%     karray_t.addCustomElement(round(element_pos(:, i), 3), 0, 2, char("el_" + string(i)));
    karray_t.addDiscElement(round(elementAll_pos(:, i), 3), 1e-5, zeros(1, 3));
%     mask = karray_t.getArrayBinaryMask(kgrid);
%     voxelPlot(double(mask))
end
mask = karray_t.getArrayBinaryMask(kgrid);
% % voxelPlot(double(mask))
imagesc(mask, [-1 1])
colormap(getColorMap);
% 
% karr2 = kWaveArray();
% karr2.addDiscElement(round(elementAll_pos(:, 1), 3), 1e-4, zeros(1, 3));
% mask2 = karr2.getArrayBinaryMask(kgrid);

end
