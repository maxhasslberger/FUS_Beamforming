function plot_results(kgrid, excitation, data, plot_title, t1_filename, plot_offset, grid_size, dx_factor, varargin)

%% Plot magnitude and phase of array elements
f = figure;
f.Position = [700 50 650 350];
subplot(2, 1, 1)
plot(abs(excitation * 1e-3))
title(plot_title);
ylabel('Amplitude (kPa)')

subplot(2, 1, 2)
plot(rad2deg(angle(excitation)))
xlabel('Element #')
ylabel('Phase (deg)')

%% Plot the pressure field 
slice_coord = 32;
dx_scan = 1e-3;

if ~isempty(varargin)
    for arg_idx = 1:2:length(varargin)
        switch varargin{arg_idx}
            case 'slice'
                slice_coord = varargin{arg_idx+1};
            case 'dx_scan'
                dx_scan = varargin{arg_idx+1};
            otherwise
                error('Unknown optional input.');
        end
    end
end

slice = round((plot_offset(2) + slice_coord) * dx_factor);

if kgrid.dim == 2
    p_data = abs(data);
else
    p_data = squeeze(abs(data(:,slice,:)));
    % Video or similar?...
end

p_sz = size(p_data);
f = figure;
f.Position = [1400 50 484 512];

% Get Ticks
plot_dy = kgrid.dy;% * 1e3;
plot_dx = kgrid.dx;% * 1e3;

plot_vecy = (plot_dy:plot_dy:grid_size(2)) / dx_scan;
plot_vecx = (plot_dx:plot_dx:grid_size(1)) / dx_scan;

plot_vecy = plot_vecy - plot_offset(3);
plot_vecx = plot_vecx - plot_offset(1);

% Include a t1w scan
if ~isempty(t1_filename)
    t1_img = niftiread(t1_filename);

    t1w_sz = [size(t1_img, 1), size(t1_img, 3)];
    if ~isequal(p_sz, t1w_sz)
        % Interpolate to adapt to grid size
        [X, Y] = meshgrid(1:t1w_sz(1), 1:t1w_sz(2));
        [Xq, Yq] = meshgrid(linspace(1, t1w_sz(1), p_sz(1)), linspace(1, t1w_sz(2), p_sz(2)));
        t1w_plot = interp2(X, Y, squeeze(double(t1_img(:, slice, :)))', Xq, Yq, "linear")';
    else
        t1w_plot = squeeze(t1_img(:, slice, :));
    end

    % Plot
    ax1 = axes;
    imagesc(ax1, plot_vecx, plot_vecy, fliplr(imrotate(t1w_plot, -90)), [50,500]);
    hold all;
    ax2 = axes;
    im2 = imagesc(ax2, plot_vecx, plot_vecy, fliplr(imrotate(p_data * 1e-3, -90)));
    im2.AlphaData = 0.5;
    linkaxes([ax1,ax2]); ax2.Visible = 'off'; ax2.XTick = []; ax2.YTick = [];
    colormap(ax1,'gray')
    colormap(ax2,'turbo')
    set([ax1,ax2],'Position',[.17 .11 .685 .815]);
    cb2 = colorbar(ax2,'Position',[.85 .11 .0275 .815]);
    xlabel(cb2, 'Pressure (kPa)');
    title(ax1,'Acoustic Pressure Amplitude')
    set(ax1, 'ydir', 'normal')
    set(ax2, 'ydir', 'normal')

else

%     plot_vecy = kgrid.y_vec(1) + kgrid.dy:kgrid.dy:kgrid.y_vec(end) + kgrid.dy;
%     plot_vecx = kgrid.x_vec(1) + kgrid.dx:kgrid.dx:kgrid.x_vec(end) + kgrid.dx;
    ax = axes;
    imagesc(ax, plot_vecx, plot_vecy, fliplr(imrotate(p_data * 1e-3, -90))); % relative to transducer 1 face (center)
    
    colormap('turbo');
    xlabel('y (mm)');
    ylabel('x (mm)');
    axis image;
    % clim([0 30000])
    title(plot_title);
    
%     c = colorbar;
%     c.Label.String = 'Pressure (kPa)';

    set(ax,'Position',[.17 .11 .685 .815]);
    cb = colorbar(ax,'Position',[.85 .11 .0275 .815]);
    xlabel(cb, 'Pressure (kPa)');
    set(ax, 'ydir', 'normal')
end

end