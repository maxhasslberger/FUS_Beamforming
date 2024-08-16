function t_mask = create_linear_array(kgrid, tr_len, x_offset, y_offset, spacing, phi)

if kgrid.dim == 2 % Create linear transducer array with approx. num_elements elements

    % Create a grid
    [Y, X] = meshgrid(1:kgrid.Ny, 1:kgrid.Nx);

    x_range = tr_len/2 * cos(phi/180*pi);
    y_range = tr_len/2 * sin(phi/180*pi);

    y_ids = [-x_range, x_range] + x_offset;
    x_ids = [-y_range, y_range] + y_offset;

    % Calculate the line equation coefficients
    A = x_ids(2) - x_ids(1);
    B = y_ids(1) - y_ids(2);
    C = y_ids(2) * x_ids(1) - y_ids(1) * x_ids(2);
    
    % Calculate the perpendicular distance from each point to the line
    distance = abs(A * X + B * Y + C) / sqrt(A^2 + B^2);

    % Find points within the threshold distance
    threshold = 0.5;
    lineMask = distance <= threshold;
    
    % Mask only the points within the bounding box and the line
    t_mask = lineMask & X >= min(y_ids) & X <= max(y_ids) & Y >= min(x_ids) & Y <= max(x_ids);
else
    error("Function not supported in 3D")
end

end