function [elementAll_pos, t_rot_per_el, active_ids, mask2el_delayFiles] = assign_arr_el_positions(elementAll_pos_orig, t_rot_per_el, active_tr_ids, n_arr_elements, n_arr_el_tot)

[elementAll_pos, el2mask_ids, mask2el_ids] = el2mask_indexing(elementAll_pos_orig, n_arr_elements);
t_rot_per_el = t_rot_per_el(:, el2mask_ids);

% Obtain transducer element ids in the right order -> t_mask
active_ids = mask2el_ids(:, active_tr_ids);
active_ids = sort(active_ids(:)); 
elementAll_pos = elementAll_pos(:, active_ids);
t_rot_per_el = t_rot_per_el(:, active_ids);

rem_el = reshape((1:n_arr_el_tot), n_arr_elements, []);
rem_el = rem_el(:, active_tr_ids);

[~, ~, mask2el_delayFiles] = el2mask_indexing(elementAll_pos_orig(:, rem_el), n_arr_elements);

end


