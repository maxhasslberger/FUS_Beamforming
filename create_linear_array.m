function t_mask = create_linear_array(kgrid, num_elements, x_offset, y_offset, spacing, rot)

if kgrid.dim == 2

    t_mask = zeros(kgrid.Nx, kgrid.Ny);

    x_range = (num_elements-1)/2 * spacing * cos(rot/180*pi);
    y_range = (num_elements-1)/2 * spacing * sin(rot/180*pi);
    
    % Determine x and y locations along the line
    xvalues = round(linspace(x_offset - x_range, x_offset + x_range, num_elements));
    yvalues = round(linspace(y_offset - y_range, y_offset + y_range, num_elements));
    
    % Replace the relevant values within the mask
    t_mask(sub2ind(size(t_mask), xvalues, yvalues)) = 1;
else
    error("Function not supported in 3D")
end

end