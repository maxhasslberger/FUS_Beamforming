function [elementAll_pos_orig, t_rot_per_el, n_arr_elements, n_arr_el_tot] = get_arr_el_positions(t_name, sparsity_name, t_pos, t_rot, t_offset_karr)

% Correct for scan offset
t_pos = t_pos + t_offset_karr;

% Init
element_pos = load(fullfile("..", "Array_Positions", t_name + ".mat")).ElementPosition'; % flat transducer array centered at [0, 0, 0] along the xy-plane

if ~strcmp(sparsity_name, "")
    sparsity_ids = load(fullfile("..", "Array_Positions", sparsity_name + ".mat")).numbers; % sparsity ids -> element_pos
    element_pos = element_pos(:, sparsity_ids);
end

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

end

