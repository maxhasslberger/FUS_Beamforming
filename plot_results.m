function plot_results(kgrid, excitation, data, t_pos, plot_title, t1_filename, t1_offset, varargin)

%% Plot magnitude and phase of array elements
figure;
subplot(2, 1, 1)
plot(abs(excitation * 1e-3))
title(plot_title);
ylabel('Amplitude (kPa)')

subplot(2, 1, 2)
plot(rad2deg(angle(excitation)))
xlabel('Element #')
ylabel('Phase (deg)')

%% Plot the pressure field 
slice_coord = 32e-3;
slice = round(t1_offset(2) + slice_coord / kgrid.dy); % Has to be overwritten

if ~isempty(varargin)
    for arg_idx = 1:2:length(varargin)
        switch varargin{arg_idx}
            case 'slice'
                slice = varargin{arg_idx+1};
            otherwise
                error('Unknown optional input.');
        end
    end
end

if kgrid.dim == 2
    p_data = abs(data);
else
    p_data = abs(data(:,:,slice));
    % Video or similar?...
end


figure;

plot_vecy = kgrid.y_vec(1) + kgrid.dy:kgrid.dy:kgrid.y_vec(end) + kgrid.dy;
plot_vecx = kgrid.x_vec(1) + kgrid.dx:kgrid.dx:kgrid.x_vec(end) + kgrid.dx;

% Include a t1w scan
if ~isempty(t1_filename)
    t1_img = niftiread(t1_filename);

    ax1 = axes;
    imagesc(ax1, imrotate(squeeze(t1_img(:, slice, :)), 90), [50,500]);
    hold all;
    ax2 = axes;
    im2 = imagesc(ax2, imrotate(p_data * 1e-3, 90));
    im2.AlphaData = 0.5;
    linkaxes([ax1,ax2]); ax2.Visible = 'off'; ax2.XTick = []; ax2.YTick = [];
    colormap(ax1,'gray')
    colormap(ax2,'turbo')
    set([ax1,ax2],'Position',[.17 .11 .685 .815]);
    cb2 = colorbar(ax2,'Position',[.85 .11 .0275 .815]);
    xlabel(cb2, 'Pressure (kPa)');
    title(ax1,'Acoustic Pressure Amplitude')
else

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

end