function plot_results(kgrid, data)

% figure;
% plot(kgrid.x_vec*1e3, amp_on_axis*1e-6, 'b.');
% set(gca, 'XLim', [0, kgrid.Nx*kgrid.dx] * 1e3);
% xlabel('Axial Position [mm]');
% ylabel('Pressure [MPa]');
% title('Axial Pressure');

if kgrid.dim == 2
    p_data = abs(data);
else
    p_data = abs(data(:,:,Nz/2+1));
    % Video or similar?...
end

% plot the pressure field 
figure;
imagesc(kgrid.x_vec, kgrid.y_vec, p_data);
colormap('turbo');
xlabel('Axial Position');
ylabel('Lateral Position');
axis image;
title('Pressure Field');

end