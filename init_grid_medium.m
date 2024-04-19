function [kgrid, medium, ppp] = init_grid_medium(f0, grid_size, varargin)

n_dim = 2;
ct_filename = [];
dx_frac = 1.0;
dx = 0;
slice_idx = 1;

if ~isempty(varargin)
    for arg_idx = 1:2:length(varargin)
        switch varargin{arg_idx}
            case 'n_dim'
                n_dim = varargin{arg_idx+1};
            case 'ct_scan'
                ct_filename = varargin{arg_idx+1};
            case 'dx_factor'
                dx_frac = varargin{arg_idx+1};
            case 'dx'
                dx = varargin{arg_idx+1};
            case 'slice_idx'
                slice_idx = varargin{arg_idx+1};
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
cfl = 0.1;

%% Define grid
% Grid size
if dx == 0
    dx = c0 / f0 / ppw;
    dx = dx * dx_frac;
end

Nx = round(grid_size(1)/dx);
Ny = round(grid_size(2)/dx);

if n_dim == 2
    kgrid = kWaveGrid(Nx, dx, Ny, dx);
    add_z = 0.0;

elseif n_dim == 3
    Nz = round(grid_size(3)/dx);

    kgrid = kWaveGrid(Nx, dx, Ny, dx, Nz, dx);
    add_z = kgrid.z_size.^2; % -> t_end

    slice_idx = 1:kgrid.Ny;
end

% Time
ppp = round(ppw / cfl); % points per temporal period
dt = 1 / (ppp * f0); % <= dx/c_max

% calculate the number of time steps to reach steady state
t_end = sqrt(kgrid.x_size.^2 + kgrid.y_size.^2 + add_z) / c0; 

Nt = round(t_end / dt);
kgrid.setTime(Nt, dt);

% t_end = 70e-6; % s
% kgrid.makeTime(medium.sound_speed, [], t_end);
% kgrid.makeTime(medium.sound_speed);

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

    skull = squeeze(skull(:, slice_idx, :));

    % assign medium properties for skull
    medium.density = rho_min + (rho_max - rho0) * ...
                    (skull - 0) / (hu_max - 0);
    medium.sound_speed = c0 + (c_max - c0) * ...
                        (medium.density - rho0) / (rho_max - rho0);
    medium.alpha_coeff = alpha_coeff_min + (alpha_coeff_max - alpha_coeff_min) * ...
                        (1 - (skull - hu_min) / (hu_max - hu_min)).^0.5;
    
    % Non-skull modeled as water
    medium.sound_speed(skull == 0) = c0;
    medium.density(skull == 0) = rho0;
    medium.alpha_coeff(skull == 0) = alpha_coeff_water;
else
    % define the homogeneous space (water)
    medium.sound_speed = c0; % * ones(Nx, Ny);
    medium.density = rho0; % * ones(Nx, Ny);
    medium.alpha_coeff = alpha_coeff_water; % dB/(MHz^y cm)
end

end
