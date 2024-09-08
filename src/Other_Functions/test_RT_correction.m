
t_name = app.ArrayElementsPositionsfilenameDropDown.Value(1:end-4);
sparsity_name = app.SparsityfilenameDropDown.Value(1:end-4);
active_tr_ids = 1:length(app.TransducerDropDown.Items);
tr_offset_karr_out = ((app.plot_offset - 1) * app.dx_scan - app.grid_size / 2)';

t_pos_orig = [-59, 68; 30, 30; 60, 60];
t_rot_orig = [0, 0; 45, -45; 180, 180];
[elementAll_pos_orig, ~, ~] = get_arr_el_positions(t_name, sparsity_name, t_pos_orig * app.dx_scan, t_rot_orig, tr_offset_karr_out);

[elementAll_pos_new, ~, n_arr_elements] = get_arr_el_positions(t_name, sparsity_name, app.t_pos * app.dx_scan, app.t_rot, tr_offset_karr_out);

err_dis = vecnorm(elementAll_pos_new - elementAll_pos_orig);
elAll_pos_new_scan = round( (elementAll_pos_new - tr_offset_karr_out) / app.kgrid.dx + app.plot_offset(:) );
elAll_pos_orig_scan = round( (elementAll_pos_orig - tr_offset_karr_out) / app.kgrid.dx + app.plot_offset(:) );
posNegDis = (( vecnorm(elAll_pos_new_scan - app.init1_id(:)) > vecnorm(elAll_pos_orig_scan - app.init1_id(:)) ) - 0.5) * 2;
err_dis = err_dis .* posNegDis;

err_time = err_dis / min(app.medium.sound_speed(:));
err_phi = err_time * 2*pi* app.CenterFreqkHzEditField.Value * 1e3;


[~, el2mask_ids_new, ~] = el2mask_indexing(elementAll_pos_new, n_arr_elements);
err_phi_sorted = err_phi(el2mask_ids_new);

[~, ~, mask2el_ids_orig] = el2mask_indexing(elementAll_pos_orig, n_arr_elements);
p_sort = app.ip.p(mask2el_ids_orig);
p_sort = p_sort(el2mask_ids_new);

app.ip.p = p_sort .* exp(1j * err_phi_sorted');
