function t_mask = create_linear_array(kgrid, tr_len, x_offset, y_offset, spacing, phi)

if kgrid.dim == 2 % Create linear transducer array with approx. num_elements elements

    t_mask = zeros(kgrid.Nx, kgrid.Ny);

    x_range = tr_len/2 * cos(phi/180*pi);
    y_range = tr_len/2 * sin(phi/180*pi);
    y_spac_factor = ((sin(phi/180*pi) >= 0.0) - 0.5) * 2;
    
    % Determine x and y locations along the line
%     xvalues = round(linspace(x_offset - x_range, x_offset + x_range, num_elements));
%     yvalues = round(linspace(y_offset - y_range, y_offset + y_range, num_elements));
    xvalues = round(x_offset - x_range):spacing:round(x_offset + x_range);
    yvalues = round(y_offset - y_range):round(y_spac_factor * spacing):round(y_offset + y_range);

    % Fit vector sizes to each other
    if length(xvalues) > length(yvalues)
        yvalues = fit_coordinates(length(xvalues), yvalues);
    elseif length(xvalues) < length(yvalues)
        xvalues = fit_coordinates(length(yvalues), xvalues);
    end
    
    disp(['xvalues: ',num2str(xvalues),'; yvalues: ',num2str(yvalues),'; size(tmask): ',num2str(size(t_mask))])

    % Replace the relevant values within the mask
    t_mask(sub2ind(size(t_mask), xvalues, yvalues)) = 1;
else
    error("Function not supported in 3D")
end

end