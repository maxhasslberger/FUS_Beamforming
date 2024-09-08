
t_name = app.ArrayElementsPositionsfilenameDropDown.Value(1:end-4);
sparsity_name = app.SparsityfilenameDropDown.Value(1:end-4);
active_tr_ids = 1:length(app.TransducerDropDown.Items);
tr_offset_karr_out = ((app.plot_offset - 1) * app.dx_scan - app.grid_size / 2)';

t_pos_orig = [-59, 68; 30, 30; 60, 60];
t_rot_orig = [0, 0; 45, -45; 180, 180];
[elementAll_pos_orig, ~, ~, ~, ~] = get_arr_el_positions(t_name, sparsity_name, t_pos_orig * app.dx_scan, t_rot_orig, tr_offset_karr_out, active_tr_ids);

t_pos_new = [-59, 68; 30, 30; 60, 60];
t_rot_new = [0, 0; 45, -45; 180, 180];
[elementAll_pos_new, ~, ~, ~, ~] = get_arr_el_positions(t_name, sparsity_name, t_pos_new * app.dx_scan, t_rot_new, tr_offset_karr_out, active_tr_ids);

err_dis = vecnorm(elementAll_pos_new - elementAll_pos_orig); %%%%%%%%%%%%%%%%%%%%%%%%%%%
err_time = err_dis / min(app.medium.sound_speed(:));
err_phi = err_time * 2*pi* app.CenterFreqkHzEditField.Value * 1e3;

app.ip.p = app.ip.p * exp(1j * err_phi);
