function [kgrid, medium, ppp] = init_grid_medium(f0, grid_size, varargin)

n_dim = 2;
ct_filename = [];
dx_factor = 1.0;
slice_idx = 1;
dx_scan = [];

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
            otherwise
                error('Unknown optional input.');
        end
    end
end

%% Constants
c0 = 1500; % m/s - water
c_max = 3100; % m/s - skull

rho0 = 1000; % kg/m^3
rho_max = 1900; % kg/m^3

alpha_coeff_water = 0.75; % dB/(MHz^y cm)
alpha_coeff_min = 4; % dB/(MHz^y cm)
alpha_coeff_max = 8.7; % dB/(MHz^y cm)

medium.alpha_power = 1.43; % Robertson et al., PMB 2017
% medium.BonA = 6; % Non-linearity

ppw = 3; % >= 2
cfl = 0.3;

%% Define grid
% Grid size
dx = c0 / f0 / ppw;
dx = dx / dx_factor;

Nx = round(grid_size(1)/dx);
Ny = round(grid_size(2)/dx);

if n_dim == 2
    kgrid = kWaveGrid(Nx, dx, Ny, dx);
    add_z = 0.0;

elseif n_dim == 3
    Nz = round(grid_size(3)/dx);

    kgrid = kWaveGrid(Nx, dx, Ny, dx, Nz, dx);
    add_z = kgrid.z_size.^2; % -> t_end

    slice_idx = 1:round(grid_size(2) / dx_scan);
end

%% Define medium
if ~isempty(ct_filename)

    hu_min = 300;
    hu_max = 2000;
    input_ct = double(niftiread(ct_filename));

    ct_max = max(input_ct(:));
    if ct_max < hu_max
        hu_max = ct_max;
    end
    
    % truncate CT HU (see Marsac et al., 2017)
    skull = input_ct;
    skull(skull < hu_min) = 0; % only use HU for skull acoustic properties
    skull(skull > hu_max) = hu_max;

    if n_dim == 2
        skull = squeeze(skull(:, slice_idx, :));
    end

    % Interpolate to adapt to grid size
    skull_sz = size(skull);
    grid_dim = size(kgrid.k);
    if ~isequal(skull_sz, grid_dim)
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
    medium.alpha_coeff = alpha_coeff_min + (alpha_coeff_max - alpha_coeff_min) * ...
                        (1 - (skull - hu_min) / (hu_max - hu_min)) .^ 0.5;
    medium.density = rho0 + (rho_max - rho0) * ...
                    (skull - 0) / (hu_max - 0);
    medium.sound_speed = c0 + (c_max - c0) * ...
                        (medium.density - rho0) / (rho_max - rho0);
    
    % Non-skull modeled as water
    medium.sound_speed(skull == 0) = c0;
    medium.density(skull == 0) = rho0;
    medium.alpha_coeff(skull == 0) = alpha_coeff_water;

    dt = cfl * dx / c_max;
else
    % define the homogeneous space (water)
    medium.sound_speed = c0; % * ones(Nx, Ny);
    medium.density = rho0; % * ones(Nx, Ny);
    medium.alpha_coeff = alpha_coeff_water; % dB/(MHz^y cm)

    dt = cfl * dx / c0;
end

% Time
ppp = round(ppw / cfl); % points per temporal period
% dt = 1 / (ppp * f0 * dx_factor);

% calculate the number of time steps to reach steady state - usually sufficient even under the presence of the skull
t_end = sqrt(kgrid.x_size.^2 + kgrid.y_size.^2 + add_z) / c0; 

Nt = round(t_end / dt);
kgrid.setTime(Nt, dt);

% t_end = 70e-6; % s
% kgrid.makeTime(medium.sound_speed, [], t_end);
% kgrid.makeTime(medium.sound_speed);

end
