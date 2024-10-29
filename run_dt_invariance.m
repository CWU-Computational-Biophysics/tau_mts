% run_dt_invariance.m
% Run a simulation (simTauMap6) at varying values of dt.

% clear the workspace
clear;

% define the save directory
save_dir = fullfile('data', 'dt_invariance');
[~, ~, ~] = mkdir(save_dir);

% iterate over three test values of 'dt'
for dt_val = [0.01, 0.001, 0.0001]
    for i = 1:20
        % load parameters
        params = genParams(ttot=5, tau_gating=true, dt=dt_val);

        % set tau to a multiple of map6
        ratio = 10;
        params.tau_on = params.map6_on * ratio;
        params.tau_off = params.tau_on;

        % print progress
        fprintf("\nStarting simulation %i @ dt=%g\n", i, params.dt)
    
        % run the simulation
        [mt_grid, mt_length, mt_state, grids, growths] = simTauMap6(params);
    
        % save the simulation
        save_name = sprintf('sim_%i_%i', dt_val*10000, i);
        save_file = fullfile(save_dir, save_name);
        save(save_file, '-nocompression', '-v7');
        fprintf("Simulation saved to '%s.mat'\n", save_file);
    end
end