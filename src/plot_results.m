function plot_results(kgrid, excitation, data, plot_title, mask2el, t1w_filename, plot_offset, grid_size, dx_factor, save_results, current_datetime, varargin)

%% Optional Inputs
slice_coord = 32;
dx_scan = 1e-3;
slice_dim = 2;
scale_factor = 1e-3;
plot_colorbar = true;

if ~isempty(varargin)
    for arg_idx = 1:2:length(varargin)
        switch varargin{arg_idx}
            case 'slice'
                slice_coord = varargin{arg_idx+1};
            case 'dx_scan'
                dx_scan = varargin{arg_idx+1};
            case 'slice_dim'
                slice_dim = varargin{arg_idx+1};
            case 'colorbar'
                plot_colorbar = varargin{arg_idx+1};
            otherwise
                error('Unknown optional input.');
        end
    end
end

%% Plot magnitude and phase of array elements
if ~isempty(excitation)
    excitation = excitation(reshape(mask2el, 1, [])); % Sort acc to transducer id
    
    f_param = figure('color','w');
    f_param.Position = [700 50 650 350];

    subplot(2, 1, 1)
    plot(abs(excitation * 1e-3), '.')
    title(plot_title);
    xlim([1 length(excitation)])
    ylabel('Amplitude (kPa)')
    
    subplot(2, 1, 2)
    plot(rad2deg(angle(excitation)), '.')
    xlim([1 length(excitation)])
    xlabel('Array Element #')
    ylabel('Phase (deg)')

    fontsize(f_param, 12,"points")
end

%% Plot the pressure field 

if kgrid.dim == 2
    p_data = abs(data);
else
    slice_p = round((plot_offset(2) + slice_coord) * dx_factor); % p space
    p_data = squeeze(abs(data(:,slice_p,:)));
end

if plot_colorbar
    cmap = turbo();
else
    cmap = hot();
end

% Get Ticks
p_sz = size(p_data);

plot_vecy = linspace(0, grid_size(2)-kgrid.dy, p_sz(2)) / dx_scan;
plot_vecx = linspace(0, grid_size(1)-kgrid.dx, p_sz(1)) / dx_scan;

plot_vecy = plot_vecy - plot_offset(3) + kgrid.dy / dx_scan;
plot_vecx = plot_vecx - plot_offset(1) + kgrid.dx / dx_scan;

f_data = figure('color','w');
f_data.Position = [1400 50 484 512];

% Include a t1w scan
if ~isempty(t1w_filename)
    t1_img = niftiread(t1w_filename);

    slice_scan = round(plot_offset(2) + slice_coord); % scan space
    t1w_sz = [size(t1_img, 1), size(t1_img, 3)];
    if ~isequal(p_sz, t1w_sz)
        % Interpolate to adapt to grid size
        [X, Y] = meshgrid(1:t1w_sz(1), 1:t1w_sz(2));
        [Xq, Yq] = meshgrid(linspace(1, t1w_sz(1), p_sz(1)), linspace(1, t1w_sz(2), p_sz(2)));
        t1w_plot = interp2(X, Y, squeeze(double(t1_img(:, slice_scan, :)))', Xq, Yq, "linear")';
    else
        t1w_plot = squeeze(t1_img(:, slice_scan, :));
    end

    % Plot
    ax1 = axes;
    imagesc(ax1, plot_vecx, plot_vecy, fliplr(imrotate(t1w_plot, -90)), [50,500]);
    hold all;
    ax2 = axes;
    im2 = imagesc(ax2, plot_vecx, plot_vecy, fliplr(imrotate(p_data * scale_factor, -90)));
    im2.AlphaData = 0.5;
    linkaxes([ax1,ax2]); ax2.Visible = 'off'; ax2.XTick = []; ax2.YTick = [];
    colormap(ax1,'gray')
    colormap(ax2,cmap)
    set([ax1,ax2],'Position',[.17 .11 .685 .815]);
    if plot_colorbar
        cb2 = colorbar(ax2,'Position',[.85 .11 .0275 .815]);
        xlabel(cb2, 'Pressure (kPa)');
    end
    title(ax1, plot_title)
    set(ax1, 'ydir', 'normal')
    set(ax2, 'ydir', 'normal')

else

    ax = axes;
    imagesc(ax, plot_vecx, plot_vecy, fliplr(imrotate(p_data * scale_factor, -90))); % relative to transducer 1 face (center)
    
    colormap(cmap);
    xlabel('x (mm)');
    ylabel('z (mm)');
    axis image;
    % clim([0 30000])
    title(plot_title);
    
%     c = colorbar;
%     c.Label.String = 'Pressure (kPa)';

    set(ax,'Position',[.17 .11 .685 .815]);
    if plot_colorbar
        cb = colorbar(ax,'Position',[.85 .11 .0275 .815]);
        xlabel(cb, 'Pressure (kPa)');
    end
    set(ax, 'ydir', 'normal')
end

fontsize(f_data, 12,"points")

if kgrid.dim == 3
    sliceViewer(double(flip(imrotate(abs(data * scale_factor), 90), 1)), 'Colormap', cmap, 'SliceNumber', slice_p, 'SliceDirection', 'Y', "Parent", figure);
    if plot_colorbar
        cb3 = colorbar;
        xlabel(cb3, 'Pressure (kPa)');
    end
    title(plot_title);
end

if save_results
    filename = fullfile("..", "Results", current_datetime + "_" + plot_title);
    export_fig(f_param, filename + "_param.pdf");
    saveas(f_data, filename + "_data.pdf");
end

end