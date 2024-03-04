function plot_results(kgrid, excitation, data, plot_title)

% Plot magnitude and phase of array elements
figure;
subplot(2, 1, 1)
plot(abs(excitation))
title(plot_title);

subplot(2, 1, 2)
plot(rad2deg(angle(excitation)))
xlabel('Element #')
ylabel('deg')

% Plot result plane
if kgrid.dim == 2
    p_data = abs(data);
else
    p_data = abs(data(:,:,Nz/2+1));
    % Video or similar?...
end

% plot the pressure field 
figure;
imagesc(p_data);
% imagesc(kgrid.x_vec, kgrid.y_vec, p_data);
colormap('turbo');
xlabel('Axial Position');
ylabel('Lateral Position');
axis image;
% clim([0 30000])
title(plot_title);
colorbar;

end