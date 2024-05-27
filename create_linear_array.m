function t_mask = create_linear_array(kgrid, num_elements, x_offset, y_offset, spacing, rot)

if kgrid.dim == 2

    t_mask = zeros(kgrid.Nx, kgrid.Ny);

    if rot == 0 || rot == 180
        start_index = round(-round(num_elements/2) * spacing + 1 + y_offset);
        t_mask(x_offset, start_index:spacing:start_index + num_elements * spacing - 1) = 1;
    elseif abs(rot) == 90
        start_index = round(-round(num_elements/2) * spacing + 1 + x_offset);
        t_mask(start_index:spacing:start_index + num_elements * spacing - 1, y_offset) = 1;
    else
        x_range = num_elements/2 * cos(rot/180*pi);
        y_range = num_elements/2 * sin(rot/180*pi);

%         num_elements = ceil(sqrt((x2 - x1).^2 + (y2 - y1).^2)) + 1;
        
        % Determine x and y locations along the line
        xvalues = round(linspace(x_offset - x_range, x_offset + x_range, num_elements));
        yvalues = round(linspace(y_offset - y_range, y_offset + y_range, num_elements));
        
        % Replace the relevant values within the mask
        t_mask(sub2ind(size(t_mask), xvalues, yvalues)) = 1;

%         tr_ids = poly2mask(round(x_offset + 0.5 * [x_range, -x_range]), round(y_offset + 0.5 * [y_range, -y_range]), kgrid.Nx, kgrid.Ny);
%         t_mask(tr_ids) = 1;
    end
else
    error("Function not supported in 3D")
end

end