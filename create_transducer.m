function karray_t = create_transducer(kgrid, t_name, t_pos, t_rot)

% Init
element_pos = load(fullfile("Array_Positions", t_name + ".mat")).ElementPosition'; % flat transducer array centered at [0, 0, 0] along the xy-plane
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
        karray_t.addDiscElement(round(elementx_pos(:, i), 3), 1e-10, ones(1, 3));
%         mask = karray_t.getArrayBinaryMask(kgrid);
%         voxelPlot(double(mask))
    end
%     mask = karray_t.getArrayBinaryMask(kgrid);
%     voxelPlot(double(mask))
end

end
