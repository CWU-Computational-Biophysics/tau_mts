% run_paper_nov26.m

% Collecting data with new parameters for publication.
% Let map_on sweep from 2.5 to 250 and map6_on = map6_off.
% tau2map6_ratio = [10, 5, 1, 0.1]

% clear the workspace
clear;

% define save directory
save_dir = fullfile('data', 'nov26_paper_run_false');
[~, ~, ~] = mkdir(save_dir);

% calculate map6 values
tau_map6_ratio = [10, 5, 1, 0.1];
default_params = genParams();
map6_vals = default_params.tau_on .* (1./tau_map6_ratio);

% iterate over values of map_on
for map6_val = map6_vals
    % repeat the simulation five times
    for i = 1:5
        % generate parameters
        params = genParams( ...
            ttot=60*10, ...
            tau_gating=false, ...
            map6_on=map6_val, ...
            map6_off=map6_val ...
            );

        % print progress
        fprintf("\nStarting simulation %i @ map6_on=%g\n", i, params.map6_on);
        fprintf("Tau Gating: %s\n", mat2str(params.tau_gating));

        % run the simulation
        [mt_grid, mt_length, mt_state, grids, growths] = simTauMap6(params);

        % save the simulation
        save_name = sprintf('sim_%g_%i.mat', params.map6_on, i);
        save_file = fullfile(save_dir, save_name);
        save(save_file, '-nocompression', '-v7');
        fprintf("Simulation saved to '%s'\n", save_file);
    end
end
