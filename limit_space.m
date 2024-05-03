function [A, b, activeA_ids] = limit_space(A, b, activeA_ids, space_limits, plot_offset, dx_factor, sound_speed)

dim = numel(size(sound_speed));
if dim == 2
    plot_offset = [plot_offset(1), plot_offset(3)];
end

space_limits = round((repmat(plot_offset', [1, 2]) + space_limits) * dx_factor);
space = zeros(size(sound_speed));

if dim == 2 % = dim
    space(space_limits(1, 1):space_limits(1, 2), space_limits(2, 1):space_limits(2, 2)) = 1;
else
    space(space_limits(1, 1):space_limits(1, 2), space_limits(2, 1):space_limits(2, 2), space_limits(3, 1):space_limits(3, 2)) = 1;
end

rem_ids = find(space);
A = A(rem_ids, :);
b = b(rem_ids);
activeA_ids = activeA_ids(rem_ids);

end