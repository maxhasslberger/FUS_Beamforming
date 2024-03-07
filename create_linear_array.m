function t_mask = create_linear_array(kgrid, num_elements, x_offset, spacing)

if kgrid.dim == 2
    
    t_mask = zeros(kgrid.Nx, kgrid.Ny);
    start_index = round(kgrid.Ny/2 - round(num_elements/2) * spacing + 1);
    t_mask(x_offset, start_index:spacing:start_index + num_elements * spacing - 1) = 1;
else
    t_mask = -1;
end