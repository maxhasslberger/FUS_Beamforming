function [karray_t, t_mask_ps, active_ids, n_arr_elements, mask2el_delayFiles] = create_transducer(kgrid, plot_offset, t_offset_karr, t_name, sparsity_name, ...
    t_pos, t_rot, active_tr_ids, el_sz)

karray_t = kWaveArray(); % real array elements including the geometry
[elementAll_pos_orig, t_rot_per_el, n_arr_elements, n_arr_el_tot] = get_arr_el_positions(t_name, sparsity_name, t_pos, t_rot, t_offset_karr);
[elementAll_pos, t_rot_per_el, active_ids, mask2el_delayFiles] = assign_arr_el_positions(elementAll_pos_orig, t_rot_per_el, active_tr_ids, n_arr_elements, n_arr_el_tot);

% Add one array element after another
t_mask_ps = zeros(size(kgrid.k));
for i = 1:length(elementAll_pos)
    new_idx = round( (elementAll_pos(:, i) - t_offset_karr) / kgrid.dx + plot_offset(:) );
    t_mask_ps(new_idx(1), new_idx(2), new_idx(3)) = 1;

    karray_t = add_array_element(karray_t, elementAll_pos(:, i), el_sz, t_rot_per_el(:, i));
end

end

