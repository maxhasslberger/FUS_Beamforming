function t_mask = create_linear_array(kgrid, num_elements, plane_offset, center_offset, spacing, rotate)

if kgrid.dim == 2
    
    t_mask = zeros(kgrid.Nx, kgrid.Ny);
    start_index = round(kgrid.Ny/2 - round(num_elements/2) * spacing + 1 + center_offset);
    t_mask(plane_offset, start_index:spacing:start_index + num_elements * spacing - 1) = 1;
    if rotate
        t_mask = rot90(t_mask);
    end
else
    error("Not supported at the moment")
end

end