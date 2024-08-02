function [kgrid, medium, sensor, sensor_mask, b_des, b_des_pl, b_mask, t_mask_ps, karray_t, only_focus_opt, ...
    active_ids, mask2el, el_per_t, t_pos, t_rot, plot_offset, point_pos, point_pos_m, grid_size, dx_factor, preplot_arg, logical_dom_ids, input_args] = ...
    init(f0, n_dim, dx_factor, varargin)

% Scan init
t1w_filename = [];
ct_filename = [];
plot_offset = [96, 127, 126] + 1; % Offset to Scan center
dx_scan = 1e-3; % m - MR Scan resolution

% Transducer and Foci init
t1_pos = [-59, 60]; % scan dims [x, z]
t2_pos = [68, 60]; % scan dims [x, z]
t_rot = [45, -45]; % deg

scan_focus_x = [-18, 22];
% scan_focus_x = [];
slice_idx_2D = 30; % Observed slice in t1w/ct scan + Ref for focus and transducer plane 
scan_focus_z = [-27, -19];
% scan_focus_z = [];
des_pressures = [300, 300]; % kPa
% des_pressures = []; % kPa

% region_labels = ["leftAmygdala", "rightAmygdala"];
% region_labels = ["leftHippocampus", "rightHippocampus"];
region_labels = [];
% des_pressures_region = [300, 300]; % kPa
des_pressures_region = []; % kPa

sidelobe_tol = 50; % percent
max_skull_pressure = 1e3; % kPa

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
    grid_size = [192, 256] * 1e-3; % m in [x, z] respectively
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
    end

    tr_offset_karr = (plot_offset * dx_scan - grid_size / 2 - kgrid.dx)'; % Compute offset for karray transducers
    grid_size = [grid_size(1), grid_size(3)]; % plane size for plots
end

dx_factor = dx_scan / round(kgrid.dx,3);

%% Segment the brain
[segment_ids] = segment_space(t1w_filename, dx_scan);
if abs(dx_factor) < 0.99 || abs(dx_factor) > 1.01
    % Interpolate to adapt to grid size
    grid_sz = size(kgrid.k);
    seg_sz = size(segment_ids);
    [uniqueStrings, ~, seg_nums] = unique(segment_ids);
    seg_nums = reshape(seg_nums, size(segment_ids)); % Ensure it has the same shape as the original 3D array
    
    % For 2D, kgrid.k has 2 dims
    [X, Y, Z] = meshgrid(1:seg_sz(1), 1:seg_sz(2), 1:seg_sz(3));
    [Xq, Yq, Zq] = meshgrid(linspace(1, seg_sz(1), grid_sz(1)), linspace(1, seg_sz(2), grid_sz(2)), linspace(1, seg_sz(3), grid_sz(3)));
    seg_nums = interp2(X, Y, Z, double(seg_nums)', Xq, Yq, Zq, "nearest")';    
    % Map back to strings
    seg_nums = round(seg_nums); % Ensure indices are integers
    segment_ids = uniqueStrings(seg_nums);
end

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

    num_elements = 75;
    % num_elements = 50;
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
%         t1_pos = [4, 55, 87]';
%         t1_rot = [45, 0, 180]'; % deg
%         t2_pos = [4, -45, 87]';
%         t2_rot = [-45, 0, 180]'; % deg
    end

    t_pos = [t1_pos, t2_pos] * 1e-3 * (1e-3 / dx_scan) + tr_offset_karr;
    % t_pos = (repmat(plot_offset', 1, size(t_pos, 2)) + t_pos) * dx_factor;
    t_rot = [t1_rot, t2_rot];
    active_tr_ids = [1, 2];

    [karray_t, t_mask_ps, active_ids, mask2el] = create_transducer(kgrid, t_name, sparsity_name, t_pos, t_rot, active_tr_ids);

    el_per_t = num_elements * ones(1, length(active_tr_ids));

%     voxelPlot(double(t_mask_ps))
end

%% Define (intracranial) Beamforming Pattern

domain_ids = segment_ids ~= "background"; % Mask entire brain
logical_dom_ids = false(numel(medium.sound_speed), 1);

cross_pixRadius = 5;

if kgrid.dim == 2
    % Mask brain slice
    slice_grid_2D = round((plot_offset(2) + slice_idx_2D) * dx_factor);
    logical_dom_ids(squeeze(domain_ids(:, slice_grid_2D, :))) = true;

    % Define targets
    if ~isempty(des_pressures)
        % Focal points - in Scan coordinate system
        point_pos_m.x = scan_focus_x;
        point_pos_m.y = scan_focus_z;
        amp_in = des_pressures' * 1e3; % Pa
    
        point_pos.x = round((plot_offset(1) + point_pos_m.x) * dx_factor);
        point_pos.y = round((plot_offset(3) + point_pos_m.y) * dx_factor);
    
        % Assign amplitude acc. to closest position
        idx = sub2ind([kgrid.Nx, kgrid.Ny], point_pos.x, point_pos.y);
        [~, order] = sort(idx);
        amp_in = amp_in(order);
    else
        point_pos_m = [];
        point_pos = [];
    end

    b_mask = zeros(kgrid.Nx, kgrid.Ny);
    amp_in_reg = des_pressures_region' * 1e3; % Pa
    point_pos.slice = slice_idx_2D;
    
    if ~only_focus_opt
        stim_regions = squeeze(segment_ids(:, slice_grid_2D, :));

        amp_vol = -1 * ones(numel(b_mask), 1);

        % Stimulate Disc pattern
        for i = 1:length(des_pressures)
            disc = makeDisc(kgrid.Nx, kgrid.Ny, point_pos.x(i), point_pos.y(i), round(0.04 * kgrid.Nx), false);
            amp_vol(logical(disc)) = amp_in(i) * ones(sum(disc(:)), 1);
            b_mask = b_mask + disc; 

            % Question ---> What is the disc pattern here and why is it added to b_mask? 
            % My understanding is that b_mask is the mask of indices of the target region 
            % and that we are using point sources, but I also noticed that the focal points seem to be stored in point_pos. 
            % Is my understanding corrent or am I off target?

            % Answer: With point_pos, we store the center points of our 
            % intended sonicated volumes (that's why point_pos.x(i) and
            % point_pos.y(i) are arguments of makeDisc). This entire part
            % is not related to the transducers and their representation as
            % point sources, but disc masks the target region around one
            % focal point, yes. At first, b_mask will be empty (line 179).
            % If you add one volume after another, you'll end up with
            % b_mask containing all the sonicated volumes. Feel free to
            % set a breakpoint in this for loop and imagesc(b_mask) after
            % each iteration. You'll see that the volumes (or areas in the
            % 2D case) will be added one after another.
        end

        % Stimulate brain region
        for j = 1:length(region_labels)
            reg_mask = stim_regions == region_labels(j);
            amp_vol(logical(reg_mask)) = amp_in_reg(j) * ones(sum(reg_mask(:)), 1);
            b_mask = b_mask + reg_mask;
        end

        b_cross = amp_vol / max(amp_vol);
        b_cross(b_cross < 0.0) = 0.0;
        amp_in = amp_vol(amp_vol >= 0);
    else

        b_cross = b_mask;
        for i = 1:length(point_pos.x)
            b_mask(point_pos.x(i), point_pos.y(i)) = 1;

            b_cross(point_pos.x(i), point_pos.y(i) - cross_pixRadius:point_pos.y(i) + cross_pixRadius) = 1;
            b_cross(point_pos.x(i) - cross_pixRadius:point_pos.x(i) + cross_pixRadius, point_pos.y(i)) = 1;
        end
    end

else
    % Mask brain
    logical_dom_ids(domain_ids) = true;
    
    % Define targets
    if ~isempty(des_pressures)
        % Focal points - in Scan coordinate system
        if isempty(t1w_filename)
            point_pos_m.x = [30, 5];
            point_pos_m.y = [10, 10];
            point_pos_m.z = [-30, 0];
        else
            point_pos_m.x = scan_focus_x;
            point_pos_m.y = slice_idx_2D * ones(1, length(scan_focus_x));
            point_pos_m.z = scan_focus_z;
        end
        amp_in = des_pressures' * 1e3; % Pa
    
        point_pos.x = round((plot_offset(1) + point_pos_m.x) * dx_factor);
        point_pos.y = round((plot_offset(2) + point_pos_m.y) * dx_factor);
        point_pos.z = round((plot_offset(3) + point_pos_m.z) * dx_factor);
    
        % Assign amplitude acc. to closest position
        idx = sub2ind([kgrid.Nx, kgrid.Ny, kgrid.Nz], point_pos.x, point_pos.y, point_pos.z);
        [~, order] = sort(idx);
        amp_in = amp_in(order);
    else
        point_pos_m = [];
        point_pos = [];
    end

    b_mask = zeros(kgrid.Nx, kgrid.Ny, kgrid.Nz);
    amp_in_reg = des_pressures_region' * 1e3; % Pa
    point_pos.slice = slice_idx_2D;

    if ~only_focus_opt

        stim_regions = segment_ids;

        amp_vol = -1 * ones(numel(b_mask), 1);

        % Stimulate Disc pattern
        for i = 1:length(des_pressures)
            ball = makeBall(kgrid.Nx, kgrid.Ny, kgrid.Nz, point_pos.x(i), point_pos.y(i), point_pos.z(i), round(0.04 * kgrid.Nx), false);
            amp_vol(logical(ball)) = amp_in(i) * ones(sum(ball(:)), 1);
            b_mask = b_mask + ball;
        end

        % Stimulate brain region
        for j = 1:length(region_labels)
            reg_mask = stim_regions == region_labels(j);
            amp_vol(logical(reg_mask)) = amp_in_reg(j) * ones(sum(reg_mask(:)), 1);
            b_mask = b_mask + reg_mask;
        end

        b_cross = amp_vol / max(amp_vol);
        b_cross(b_cross < 0.0) = 0.0;
        amp_in = amp_vol(amp_vol >= 0);
    else
    
        b_cross = b_mask;
        for i = 1:length(des_pressures)
            b_mask(point_pos.x(i), point_pos.y(i), point_pos.z(i)) = 1;

            b_cross(point_pos.x(i), point_pos.y(i), point_pos.z(i) - cross_pixRadius:point_pos.z(i) + cross_pixRadius) = 1;
            b_cross(point_pos.x(i), :, point_pos.z(i)) = 1;
            b_cross(point_pos.x(i) - cross_pixRadius:point_pos.x(i) + cross_pixRadius, point_pos.y(i), point_pos.z(i)) = 1;
        end
    end

end
b_mask = logical(b_mask);

% Create preview plot
preplot_arg = reshape(b_cross, size(b_mask));
preplot_arg(logical(t_mask_ps)) = max(b_cross(:));

if ~isscalar(medium.sound_speed)
    skull_arg = medium.sound_speed / max(medium.sound_speed(:));
    skull_arg = skull_arg - min(skull_arg(:));
    preplot_arg = preplot_arg + skull_arg * max(b_cross(:));
end

plot_results(kgrid, [], preplot_arg, 'Plot Preview', [], t1w_filename, plot_offset, grid_size, dx_factor, false, [], 'slice', point_pos.slice, 'colorbar', false, ...
    'cmap', hot());
preplot_arg(logical(t_mask_ps)) = 0.0; % Do not show transducers in second pre-plot

% Create desired signal
phase = zeros(length(amp_in), 1); % Zero phase for entire observation plane

b_des = amp_in .* exp(1j*phase); % only observed elements

b_max = max(abs(b_des));
b_des_pl = sidelobe_tol/100 * b_max * ones(numel(b_mask), 1); % Entire plane max amp
b_des_pl(b_mask) = b_des; % Target amp
b_des_pl(logical(reshape(medium.sound_speed > min(medium.sound_speed(:)), [], 1)) & ~logical_dom_ids) = max_skull_pressure * 1e3; % Skull max amp

% set simulation input options
input_args = {'PMLSize', 10, 'PMLInside', true, 'PlotPML', true, 'DisplayMask', b_mask + t_mask_ps, 'RecordMovie', false};

end