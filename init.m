function [kgrid, medium, sensor, sensor_mask, b_des, b_des_pl, b_mask, t_mask_ps, karray_t, only_focus_opt, space_limits, ...
    active_ids, mask2el, el_per_t, t_pos, t_rot, plot_offset, point_pos, point_pos_m, grid_size, dx_factor, input_args] = ...
    init(f0, n_dim, dx_factor, varargin)

% Scan init
t1w_filename = [];
ct_filename = [];
plot_offset = [96, 127, 126] + 1; % Offset to Scan center
dx_scan = 1e-3; % m - Scan resolution

% Transducer and Foci init
t1_pos = [-40, 79]; % scan dims
t2_pos = [47, 79]; % scan dims
t_rot = [30, -30]; % deg

scan_focus_x = [-18, 18];
slice_idx_2D = 30; % Observed slice in t1w/ct scan - Serves as focus plane as well
scan_focus_z = [-27, -18];
des_pressures = [300, 300]; % kPa

sidelobe_tol = 10; % percent

% Simulation config
only_focus_opt = true; % Optimize only focal spots or entire grid
use_greens_fctn = true; % Use Green's function to obtain propagation matrix A (assuming point sources and a lossless homogeneous medium)

if ~isempty(varargin)
    for arg_idx = 1:2:length(varargin)
        switch varargin{arg_idx}
            case 'sidelobe_tol'
                sidelobe_tol = varargin{arg_idx+1};
            case 't1_scan'
                t1w_filename = varargin{arg_idx+1};
            case 'ct_scan'
                ct_filename = varargin{arg_idx+1};
            case 'slice_idx'
                slice_idx_2D = varargin{arg_idx+1};
            case 'dx_scan'
                dx_scan = varargin{arg_idx+1};
            case 'only_focus_opt'
                only_focus_opt = varargin{arg_idx+1};
            case 'use_greens_fctn'
                use_greens_fctn = varargin{arg_idx+1};
            otherwise
                error('Unknown optional input.');
        end
    end
end

if n_dim == 2
    grid_size = [192, 256] * 1e-3; % m in [x, y] respectively
else
    if isempty(t1w_filename)
        grid_size = [140, 100, 140] * 1e-3; % m in [x, y, z] respectively
%         grid_size = [192, 256, 256] * 1e-3; % m in [x, y, z] respectively
    else
        grid_size = [192, 256, 256] * 1e-3; % m in [x, y, z] respectively
    end
end

[kgrid, medium, ppp] = init_grid_medium(f0, grid_size, 'n_dim', n_dim, 'dx_factor', dx_factor, 'ct_scan', ct_filename, ...
    'slice_idx', round(plot_offset(2) + slice_idx_2D), 'dx_scan', dx_scan);
[sensor, sensor_mask] = init_sensor(kgrid, ppp);

if n_dim == 3
    if isempty(t1w_filename)
        plot_offset = grid_size / dx_scan / 2; % Offset to center
%         dx_scan = kgrid.dx;
    end

    tr_offset_3D = (plot_offset * dx_scan - grid_size / 2 - kgrid.dx)';
    grid_size = [grid_size(1), grid_size(3)];
end

dx_factor = dx_scan / kgrid.dx;

%% Define Transducer Geometry

if kgrid.dim == 2
    % Transducer coordinates and alignment- in Scan coordinate system
%     t1_pos = [-70, -2]'; % scan dims
%     t2_pos = [0, 83]'; % scan dims
%     t3_pos = [77, -2]'; % scan dims
%     t_rot = [90, 0, 90];

    t1_pos = t1_pos'; % scan dims
    t2_pos = t2_pos'; % scan dims
    t_rot = t_rot';

    t_pos = [t1_pos, t2_pos];

    num_elements = 40;
    spacing = ceil(dx_factor);
    n_trs = length(t_rot);

    t_mask_ps = false(kgrid.Nx, kgrid.Ny);
    el_per_t = zeros(1, n_trs);
    t_ids = [];
    for i = 1:n_trs
        el_offset = round((plot_offset(1) + t_pos(1, i)) * dx_factor); % grid points
        shift = round((plot_offset(3) + t_pos(2, i)) * dx_factor); % tangential shift in grid points
    
        new_arr = create_linear_array(kgrid, num_elements, el_offset, shift, spacing, t_rot(i));

        el_per_t(i) = sum(new_arr(:));
        t_ids = [t_ids; find(new_arr)];
        t_mask_ps = t_mask_ps | logical(new_arr);
    end
    
    [~, el2mask_ids] = sort(t_ids);
    [~, mask2el] = sort(el2mask_ids);

    karray_t = [];
    active_ids = [];

%     imagesc(t_mask_ps, [-1 1])
%     colormap(getColorMap);
else
    % Planar Array
    if use_greens_fctn
        t_name = "std";
    else
        t_name = "std_orig";
    end
%     t_name = "std_orig";

    sparsity_name = "sparsity_ids";
    num_elements = 128;

    % Transducer coordinates and alignment- in Scan coordinate system
    if isempty(t1w_filename)
        t1_pos = [30, 0, 65]';
        t1_rot = [0, 0, 180]'; % deg
        t2_pos = [-65, 0, -30]';
        t2_rot = [180, 90, 0]'; % deg
    else
        t1_pos = [t1_pos(1), slice_idx_2D, t1_pos(2)]';
        t1_rot = [0, t_rot(1), 180]'; % deg
        t2_pos = [t2_pos(1), slice_idx_2D, t2_pos(2)]';
        t2_rot = [0, t_rot(2), 180]'; % deg
    end

    t_pos = [t1_pos, t2_pos] * 1e-3 * (1e-3 / dx_scan) + tr_offset_3D;
    % t_pos = (repmat(plot_offset', 1, size(t_pos, 2)) + t_pos) * dx_factor;
    t_rot = [t1_rot, t2_rot];
    active_tr_ids = [1, 2];

    [karray_t, t_mask_ps, active_ids, mask2el] = create_transducer(kgrid, t_name, sparsity_name, t_pos, t_rot, active_tr_ids);

    el_per_t = num_elements * ones(1, length(active_tr_ids));

%     voxelPlot(double(t_mask_ps))
end

%% Define (intracranial) Beamforming Pattern
cross_pixRadius = 5;

if kgrid.dim == 2

    % Focal points - in Scan coordinate system
    point_pos_m.x = scan_focus_x;
    point_pos.slice = slice_idx_2D;
    point_pos_m.y = scan_focus_z;
    amp_in = des_pressures' * 1e3; % Pa

    point_pos.x = round((plot_offset(1) + point_pos_m.x) * dx_factor);
    point_pos.y = round((plot_offset(3) + point_pos_m.y) * dx_factor);

    % Assign amplitude acc. to closest position
    idx = sub2ind([kgrid.Nx, kgrid.Ny], point_pos.x, point_pos.y);
    [~, order] = sort(idx);
    amp_in = amp_in(order);

    b_mask = zeros(kgrid.Nx, kgrid.Ny);
    b_cross = b_mask;
    
    if ~only_focus_opt

        % Stimulate Disc pattern
        for i = 1:length(point_pos.x)
            b_mask = b_mask + makeDisc(kgrid.Nx, kgrid.Ny, point_pos.x(i), point_pos.y(i), round(0.025 * kgrid.Nx), false);
            amp_in = amp_in(i) * ones(sum(b_mask(:)), 1);
        end

        b_cross = b_mask;
        space_limits = [-67, 71; -74, 86];
    else

        for i = 1:length(point_pos.x)
            b_mask(point_pos.x(i), point_pos.y(i)) = 1;
            b_cross(point_pos.x(i), point_pos.y(i) - cross_pixRadius:point_pos.y(i) + cross_pixRadius) = 1;
            b_cross(point_pos.x(i) - cross_pixRadius:point_pos.x(i) + cross_pixRadius, point_pos.y(i)) = 1;
        end

        space_limits = [];
    end

else
    only_focus_opt = true;
    space_limits = [];

    % Focal points - in Scan coordinate system
    if isempty(t1w_filename)
        point_pos_m.x = [30, 5];
        point_pos_m.y = [10, 10];
        point_pos_m.z = [-30, 0];
        amp_in = des_pressures' * 1e3; % Pa
    else
        point_pos_m.x = scan_focus_x;
        point_pos_m.y = [slice_idx_2D, slice_idx_2D];
        point_pos_m.z = scan_focus_z;
        amp_in = des_pressures' * 1e3; % Pa
    end

    point_pos.x = round((plot_offset(1) + point_pos_m.x) * dx_factor);
    point_pos.y = round((plot_offset(2) + point_pos_m.y) * dx_factor);
    point_pos.z = round((plot_offset(3) + point_pos_m.z) * dx_factor);

    % Assign amplitude acc. to closest position
    idx = sub2ind([kgrid.Nx, kgrid.Ny, kgrid.Nz], point_pos.x, point_pos.y, point_pos.z);
    [~, order] = sort(idx);
    amp_in = amp_in(order);

    b_mask = zeros(kgrid.Nx, kgrid.Ny, kgrid.Nz);

    for point = 1:length(point_pos.x)
        b_mask(point_pos.x(point), point_pos.y(point), point_pos.z(point)) = 1;
    end

    point_pos.slice = point_pos_m.y(1);

    b_cross = b_mask;
    for i = 1:length(point_pos.x)
        b_cross(point_pos.x(i), point_pos.y(i), point_pos.z(i) - cross_pixRadius:point_pos.z(i) + cross_pixRadius) = 1;
        b_cross(point_pos.x(i), :, point_pos.z(i)) = 1;
        b_cross(point_pos.x(i) - cross_pixRadius:point_pos.x(i) + cross_pixRadius, point_pos.y(i), point_pos.z(i)) = 1;
    end

end

% Create preview plot
if ~isscalar(medium.sound_speed)
    plot_arg = medium.sound_speed / max(medium.sound_speed(:));
    plot_arg = plot_arg - min(plot_arg(:));
else
    plot_arg = zeros(size(t_mask_ps));
end

plot_arg(logical(t_mask_ps)) = 0.5;
plot_arg(logical(b_cross)) = 1.0;

plot_results(kgrid, [], plot_arg, 'Plot Preview', [], t1w_filename, plot_offset, grid_size, dx_factor, false, [], 'slice', point_pos.slice)

% Create desired signal
phase = zeros(length(amp_in), 1); % Zero phase for entire observation plane

b_des = amp_in .* exp(1j*phase); % only observed elements

b_max = max(abs(b_des));
b_des_pl = sidelobe_tol/100 * b_max * ones(kgrid.Nx * kgrid.Ny, 1); % entire plane
b_des_pl(logical(b_mask)) = b_des;

% set simulation input options
input_args = {'PMLSize', 10, 'PMLInside', true, 'PlotPML', true, 'DisplayMask', b_mask + t_mask_ps, 'RecordMovie', false};

end