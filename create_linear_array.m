function t_mask = create_linear_array(kgrid, num_elements, x_offset, y_offset, spacing, rotate90)
% TODO: Use karray class (-> arbitrary rotation)
if kgrid.dim == 2
    
    t_mask = zeros(kgrid.Nx, kgrid.Ny);
    if ~rotate90
        start_index = round(-round(num_elements/2) * spacing + 1 + y_offset);
        t_mask(x_offset, start_index:spacing:start_index + num_elements * spacing - 1) = 1;
    else
        start_index = round(-round(num_elements/2) * spacing + 1 + x_offset);
        t_mask(start_index:spacing:start_index + num_elements * spacing - 1, y_offset) = 1;
    end
else
    error("Not supported at the moment")
end

end