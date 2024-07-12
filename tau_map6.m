% tau_map6.m
% Tau/Map6 project coarse grain model

% define user-input simulation parameters
% total simulation time [second]
ttot = 50;
% T0/M0 ratio [#]
tau_map_ratio = 50;
% force leading edge tau binding [bool]
force_tau = false;

% grid parameters
% initial number of spatial grid points [#]
grid_points = 20;
% funamental unit of length [micrometer]
dx = 0.036;

% time parameters
% time scale [second]
dt = 0.001;
% total number of time steps [#]
steps = ttot/dt;

% initialize simple arrays
% time array [second]
t = 0:dt:steps;
% position array [micrometer]
x = 0:grid_points;

% initialize the MT grid
% 0 = empty, 1 = tau, 2 = map6
nonexist = 0;
empty = 1;
tau = 2;
map6 = 3;
mt_grid = zeros(steps, grid_points);