% run.m
% Run a simulation (simTauMap6) and save the results.

% clear the workspace
clear;

% define save file name
save_name = 'simulation';

% initialize simulation parameters
params = genParams(ttot=10, tau_gating=true);

% set tau to a multiple of map6
ratio = 10;
params.tau_on = params.map6_on * ratio;
params.tau_off = params.tau_on;


% run the simulation
[mt_grid, mt_length, mt_state, grids, growths] = simTauMap6(params);


% save the simulation
save_dir = fullfile('data', 'tests');
[~, ~, ~] = mkdir(save_dir);
save_file = fullfile(save_dir, save_name);
save(save_file, '-nocompression', '-v7');
fprintf("Simulation saved to '%s.mat'\n", save_file)