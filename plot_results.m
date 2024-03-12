function plot_results(kgrid, excitation, data, plot_title)

% Plot magnitude and phase of array elements
figure;
subplot(2, 1, 1)
plot(abs(excitation * 1e-3))
title(plot_title);
ylabel('Amplitude (kPa)')

subplot(2, 1, 2)
plot(rad2deg(angle(excitation)))
xlabel('Element #')
ylabel('Phase (deg)')

% Plot result plane
if kgrid.dim == 2
    p_data = abs(data);
else
    p_data = abs(data(:,:,Nz/2+1));
    % Video or similar?...
end

% plot the pressure field 
figure;
imagesc(p_data * 1e-3);
% imagesc(kgrid.x_vec, kgrid.y_vec, p_data);
colormap('turbo');
xlabel('x Position (grid points)');
ylabel('y Position (grid points)');
axis image;
% clim([0 30000])
title(plot_title);
c = colorbar;
c.Label.String = 'Pressure (kPa)';

end