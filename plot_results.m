function plot_results(kgrid, excitation, data, t_pos, plot_title, varargin)

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
    z_coord = round(kgrid.Nz/2);

    if ~isempty(varargin)
        for arg_idx = 1:2:length(varargin)
            switch varargin{arg_idx}
                case 'z_coord'
                    z_coord = varargin{arg_idx+1};
                otherwise
                    error('Unknown optional input.');
            end
        end
    end

    p_data = abs(data(:,:,z_coord));
    % Video or similar?...
end

% plot the pressure field 
figure;

plot_vecy = kgrid.y_vec(1):kgrid.dy:kgrid.y_vec(end);
plot_vecx = kgrid.x_vec(1):kgrid.dx:kgrid.x_vec(end);

imagesc((plot_vecy- t_pos(2, 1)) * 1e3, (plot_vecx - t_pos(1, 1)) * 1e3, p_data * 1e-3); % relative to transducer 1 face (center)

colormap('turbo');
xlabel('y (mm)');
ylabel('x (mm)');
axis image;
% clim([0 30000])
title(plot_title);

c = colorbar;
c.Label.String = 'Pressure (kPa)';

end