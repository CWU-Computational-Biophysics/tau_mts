% run_dt_invariance.m
% Run a simulation (simTauMap6) at varying values of dt.

% clear the workspace
clear;

% define the save directory
save_dir = fullfile('data', 'tm_sweep');
[~, ~, ~] = mkdir(save_dir);

% iterate over four ratios
for ratio_val = [1, 5, 10, 50]
    for i = 1:5
        % load parameters
        params = genParams(ttot=20, tau_gating=true);

        % set tau to a multiple of map6
        ratio = ratio_val;
        params.tau_on = params.map6_on * ratio;
        params.tau_off = params.tau_on;

        % print progress
        fprintf("\nStarting simulation %i @ tau_on=%g\n", i, params.tau_on)

        % run the simulation
        [mt_grid, mt_length, mt_state, grids, growths] = simTauMap6(params);

        % save the simulation
        save_name = sprintf('sim_%i_%i', ratio_val, i);
        save_file = fullfile(save_dir, save_name);
        save(save_file, '-nocompression', '-v7');
        fprintf("Simulation saved to '%s.mat'\n", save_file);
    end
end
