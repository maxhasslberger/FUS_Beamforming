function [karray_t, t_mask_ps, active_ids, mask2el_delayFiles] = create_transducer(kgrid, t_name, sparsity_name, t_pos, t_rot, active_tr_ids)

% Init
element_pos = load(fullfile("..", "Array_Positions", t_name + ".mat")).ElementPosition'; % flat transducer array centered at [0, 0, 0] along the xy-plane
sparsity_ids = load(fullfile("..", "Array_Positions", sparsity_name + ".mat")).numbers; % sparsity ids -> element_pos
el_sz = [3, 3] * 1e-3; % dimensions of one array element

element_pos = element_pos(:, sparsity_ids);

karray_t = kWaveArray(); % real array elements including the geometry
karray_t_ps = kWaveArray(); % approximation with point sources for optimization
n_arr_elements = size(element_pos, 2);
n_arr_el_tot = n_arr_elements * size(t_pos, 2);

% Obtain std configuration with all transducers
elementAll_pos_orig = nan(3, n_arr_el_tot);
t_rot_per_el = nan(3, n_arr_el_tot);
for tx = 1:size(t_pos, 2) % for each transducer

    tx_pos = t_pos(:, tx);
    tx_rot = t_rot(:, tx);

    % Shift and align transducers 
    tx_transf = getAffineMatrix(tx_pos, tx_rot);
    
    elementx_pos = tx_transf(1:3, 1:3) * element_pos;
    elementAll_pos_orig(:, (tx - 1) * n_arr_elements + 1 : tx * n_arr_elements) = elementx_pos + repmat(tx_pos, 1, size(element_pos, 2));

    t_rot_per_el(:, (tx - 1) * n_arr_elements + 1 : tx * n_arr_elements) = repmat(tx_rot, 1, size(element_pos, 2));
end

[elementAll_pos, el2mask_ids, mask2el_ids] = el2mask_indexing(elementAll_pos_orig, n_arr_elements);
t_rot_per_el = t_rot_per_el(:, el2mask_ids);

% Obtain transducer element ids in the right order -> t_mask
active_ids = mask2el_ids(:, active_tr_ids);
active_ids= sort(active_ids(:)); 
elementAll_pos = elementAll_pos(:, active_ids);
t_rot_per_el = t_rot_per_el(:, active_ids);

rem_el = reshape((1:n_arr_el_tot), n_arr_elements, []);
rem_el = rem_el(:, active_tr_ids);

[~, ~, mask2el_delayFiles] = el2mask_indexing(elementAll_pos_orig(:, rem_el), n_arr_elements); 

% Add one array element after another
for i = 1:length(elementAll_pos)

    karray_t_ps.addDiscElement(elementAll_pos(:, i), 1e-10, zeros(1, 3));
    karray_t.addRectElement(elementAll_pos(:, i), el_sz(1), el_sz(2), t_rot_per_el(:, i));

%     mask = karray_t_ps.getArrayBinaryMask(kgrid);
%     voxelPlot(double(mask))
end

t_mask_ps = karray_t_ps.getArrayBinaryMask(kgrid);

% voxelPlot(double(t_mask_ps))

% t_mask = karray_t.getArrayBinaryMask(kgrid);
% voxelPlot(double(t_mask))

% figure
% sliceViewer(karray_t.getArrayGridWeights(kgrid))

end


function [elementAll_pos, el2mask_ids, mask2el_ids] = el2mask_indexing(elementAll_pos_orig, n_arr_elements)

elementAll_pos = flip(elementAll_pos_orig, 1);

[elementAll_pos, el2mask_ids] = sortrows(elementAll_pos'); % refer element indices to mask -> right order in getDistributedSourceSignal

elementAll_pos = elementAll_pos';
elementAll_pos = flip(elementAll_pos, 1);

[~, mask2el_ids] = sort(el2mask_ids); % Refer back to original order
mask2el_ids = reshape(mask2el_ids, n_arr_elements, []); % Original order per transducer assuming each transducer has the same amount of elements

end
