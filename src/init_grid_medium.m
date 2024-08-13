function [kgrid, medium, grid_size, ppp] = init_grid_medium(f0, grid_size, varargin)

n_dim = 2;
ct_filename = [];
dx_factor = 1.0;
slice_idx = 1;
dx_scan = 1e-3;
const = [];

if ~isempty(varargin)
    for arg_idx = 1:2:length(varargin)
        switch varargin{arg_idx}
            case 'n_dim'
                n_dim = varargin{arg_idx+1};
            case 'ct_scan'
                ct_filename = varargin{arg_idx+1};
            case 'dx_factor'
                dx_factor = varargin{arg_idx+1};
            case 'slice_idx'
                slice_idx = varargin{arg_idx+1};
            case 'dx_scan'
                dx_scan = varargin{arg_idx+1};
            case 'constants'
                const = varargin{arg_idx+1};
            otherwise
                error('Unknown optional input.');
        end
    end
end

%% Constants
if isempty(const)
    const.c0 = 1500; % m/s - water
    const.c_max = 3100; % m/s - skull
    
    const.rho0 = 1000; % kg/m^3
    const.rho_max = 1900; % kg/m^3
    
    const.alpha_coeff_water = 0.75; % dB/(MHz^y cm)
    const.alpha_coeff_min = 4; % dB/(MHz^y cm)
    const.alpha_coeff_max = 8.7; % dB/(MHz^y cm)
    
    const.alpha_power = 1.43; % Robertson et al., PMB 2017
    % medium.BonA = 6; % Non-linearity
    
    const.hu_min = 300; % Hounsfield units
    const.hu_max = 2000;
    
    const.ppw = 3; % >= 2
    const.cfl = 0.3;
end

medium.alpha_power = const.alpha_power;

%% Define grid
% Grid size
dx = const.c0 / f0 / const.ppw;
dx = dx / dx_factor;

if ~isempty(ct_filename)
    input_ct = double(niftiread(ct_filename));
    skull = input_ct;
    skull_sz = size(skull);

    skull_res_factor = dx_scan / dx;
    Nx = round(skull_sz(1) * skull_res_factor);
    Ny = round(skull_sz(2) * skull_res_factor);
    Nz = round(skull_sz(3) * skull_res_factor);

    grid_size = [Nx, Ny, Nz] * dx;
else
    Nx = round(grid_size(1)/dx);
    Ny = round(grid_size(2)/dx);
end

if n_dim == 2
    kgrid = kWaveGrid(Nx, dx, Ny, dx);
    add_z = 0.0;

elseif n_dim == 3
    if ~exist("Nz", "var")
        Nz = round(grid_size(3)/dx);
    end

    kgrid = kWaveGrid(Nx, dx, Ny, dx, Nz, dx);
    add_z = kgrid.z_size.^2; % -> t_end

    slice_idx = 1:round(grid_size(2) / dx_scan);
end

%% Define medium
if ~isempty(ct_filename)

    ct_max = max(input_ct(:));
    if ct_max < const.hu_max
        const.hu_max = ct_max;
    end
    
    % truncate CT HU (see Marsac et al., 2017)
    skull(skull < const.hu_min) = 0; % only use HU for skull acoustic properties
    skull(skull > const.hu_max) = const.hu_max;

    if n_dim == 2
        skull = squeeze(skull(:, slice_idx, :));
    end

    % Interpolate to adapt to grid size
    grid_dim = size(kgrid.k);
%     if ~isequal(skull_sz, grid_dim)
    if ~isequal(dx, dx_scan)
        if n_dim == 2
            [X, Y] = meshgrid(1:skull_sz(1), 1:skull_sz(2));
            [Xq, Yq] = meshgrid(linspace(1, skull_sz(1), grid_dim(1)), linspace(1, skull_sz(2), grid_dim(2)));
            skull = interp2(X, Y, skull', Xq, Yq, "linear")';
        else
            [X, Y, Z] = meshgrid(1:skull_sz(1), 1:skull_sz(2), 1:skull_sz(3));
            [Xq, Yq, Zq] = meshgrid(linspace(1, skull_sz(1), grid_dim(1)), linspace(1, skull_sz(2), grid_dim(2)), linspace(1, skull_sz(3), grid_dim(3)));
            skull = permute(interp3(X, Y, Z, permute(skull, [2 1 3]), Xq, Yq, Zq, "linear"), [2 1 3]);
        end
    end

    % assign medium properties for skull
    medium.alpha_coeff = const.alpha_coeff_min + (const.alpha_coeff_max - const.alpha_coeff_min) * (1 - (skull - const.hu_min) / (const.hu_max - const.hu_min)) .^ 0.5;
    medium.density = const.rho0 + (const.rho_max - const.rho0) * (skull - 0) / (const.hu_max - 0);
    medium.sound_speed = const.c0 + (const.c_max - const.c0) * (medium.density - const.rho0) / (const.rho_max - const.rho0);
    
    % Non-skull modeled as water
    medium.sound_speed(skull == 0) = const.c0;
    medium.density(skull == 0) = const.rho0;
    medium.alpha_coeff(skull == 0) = const.alpha_coeff_water;

    dt = const.cfl * dx / const.c_max;
else
    % define the homogeneous space (water)
    medium.sound_speed = const.c0; % * ones(Nx, Ny);
    medium.density = const.rho0; % * ones(Nx, Ny);
    medium.alpha_coeff = const.alpha_coeff_water; % dB/(MHz^y cm)

    dt = const.cfl * dx / const.c0;
end

% Time
ppp = round(1 / (dt * f0)); % points per temporal period
% ppp = round(ppw / cfl); % points per temporal period
% dt = 1 / (ppp * f0);
% dt_stability_lim = checkStability(kgrid, medium);
% if dt_stability_lim ~= Inf
%     dt = dt_stability_lim;
% end

% calculate the number of time steps to reach steady state - usually sufficient even under the presence of the skull
t_end = sqrt(kgrid.x_size.^2 + kgrid.y_size.^2 + add_z) / const.c0; 

Nt = round(t_end / dt);
kgrid.setTime(Nt, dt);

end
