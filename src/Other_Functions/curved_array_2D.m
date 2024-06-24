% Curvature
grid_size = [kgrid.Nx, kgrid.Ny];
radius = kgrid.Nx; % grid points -> curvature
diameter = round(0.6 * kgrid.Ny + 1 - mod(kgrid.Ny, 2)); % grid points -> one end to the other
focus_pos = round(grid_size / 2);

t_mask = makeArc(grid_size, [el1_offset, round(kgrid.Ny / 2 + shift)], radius, diameter, focus_pos);
t_mask = t_mask + makeArc(grid_size, [round(kgrid.Nx / 2 + shift), el2_offset], radius, diameter, focus_pos); % Second (orthogonal) curved array
% t_mask = t_mask + makeArc(grid_size, [round(kgrid.Nx - el1_offset), round(kgrid.Ny / 2 + shift)], radius, diameter, focus_pos); % Second (antiparallel) curved array
