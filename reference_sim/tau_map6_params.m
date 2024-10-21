% tau_map6_params.m
% Tau/Map6 simulation user-input parameters

% define user-input simulation parameters
% total simulation time [second]
ttot = 50;
% T0/M0 ratio [#]
tm_ratio = 10;
% force leading edge tau binding [bool]
t_force = true;
% simulation save name
datetime_str = string(datetime("now", Format="uuuuMMMdd-HHmmss"));
sim_name = sprintf("taumap6_" + datetime_str);
export_dir = fullfile("data/tests");

% grid parameters
% initial number of spatial grid points [#] (M=20)
grid_points_init = 20;
% funamental unit of length [micrometer] (dx=0.036)
dx = 0.036;
% time scale [second] (dt=0.001)
dt = 0.01;

% tunable parameters
% base binding rate for map6 [1/second] (M0=0.1)
m_on = 0.1;
% unbinding rate for map6 [1/second] (Moff=0.1)
m_off = m_on;
% cooperativity for tau [ ] (alphaT=0.0)
alpha_t = 0.0;
% cooperativity for map6 [ ] (alphaM=0.0)
alpha_m = 0.0;
% rescue frequency [1/second] (fmp=500)
f_res = 50;
% catastrophe frequency [1/second] (fpm=300)
f_cat = 30;
