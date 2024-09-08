function [karray_t, t_mask_ps, active_ids, num_elements, mask2el_delayFiles] = create_transducer(kgrid, plot_offset, t_offset_karr, t_name, sparsity_name, ...
    t_pos, t_rot, active_tr_ids, el_sz)

karray_t = kWaveArray(); % real array elements including the geometry
[elementAll_pos, t_rot_per_el, active_ids, num_elements, mask2el_delayFiles] = get_arr_el_positions(t_name, sparsity_name, t_pos, t_rot, t_offset_karr, active_tr_ids); 

% Add one array element after another
t_mask_ps = zeros(size(kgrid.k));
for i = 1:length(elementAll_pos)
    new_idx = round( (elementAll_pos(:, i) - t_offset_karr) / kgrid.dx + plot_offset(:) );
    t_mask_ps(new_idx(1), new_idx(2), new_idx(3)) = 1;

    karray_t = add_array_element(karray_t, elementAll_pos(:, i), el_sz, t_rot_per_el(:, i));
end

end

