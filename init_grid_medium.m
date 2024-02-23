function [kgrid, medium, ppp] = init_grid_medium(f0, varargin)

n_dim = 2;

if ~isempty(varargin)
    for arg_idx = 1:2:length(varargin)
        switch varargin{arg_idx}
            case 'n_dim'
                n_dim = varargin{arg_idx+1};
            case 'ct_scan'
                ct_scan = varargin{arg_idx+1};
            otherwise
                error('Unknown optional input.');
        end
    end
end

%% Constants
c0 = 1500; % m/s - water
rho0 = 1000; % kg/m^3

ppw = 3; % >= 2
cfl = 0.3;

%% Define grid

% Grid size
x_size = 100e-3; % m
y_size = x_size;

dx = c0 / f0 / ppw;

Nx = round(x_size/dx);
Ny = round(y_size/dx);
N = [Nx, Ny];

if n_dim == 2
    kgrid = kWaveGrid(Nx, dx, Ny, dx);
    add_z = 0.0;

elseif n_dim == 3
    z_size = x_size;
    Nz = round(z_size/dx);
    N = [N, Nz];

    kgrid = kWaveGrid(Nx, dx, Ny, dx, Nz, dx);
    add_z = kgrid.z_size.^2;
end

% Time
ppp = round(ppw / cfl); % points per temporal period
dt = 1 / (ppp * f0);

% calculate the number of time steps to reach steady state
t_end = sqrt(kgrid.x_size.^2 + kgrid.y_size.^2 + add_z) / c0; 

Nt = round(t_end / dt);
kgrid.setTime(Nt, dt);

% t_end = 70e-6; % s
% kgrid.makeTime(medium.sound_speed, [], t_end);
% kgrid.makeTime(medium.sound_speed);

%% Define medium

% define the properties of the propagation medium
medium.sound_speed = c0; % * ones(Nx, Ny);
medium.density = rho0; % * ones(Nx, Ny);

% medium.alpha_coeff = 0.75; % dB/(MHz^y cm)
% medium.alpha_power = 1.43; % Robertson et al., PMB 2017
% medium.BonA = 6;

end
