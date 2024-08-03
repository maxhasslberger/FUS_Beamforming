function short_arr = fit_coordinates(len_long, short_arr)
% Copy array and fill up remaining indices
short_mat = repmat(short_arr, floor(len_long / length(short_arr)), 1);
rem_quant = mod(len_long, length(short_arr));
add_short = [short_arr(1:rem_quant), nan(1, length(short_arr) - rem_quant)];

% Transform to 1D
short_mat = [short_mat; add_short];
short_arr = short_mat(:)';
short_arr(isnan(short_arr)) = [];

end
