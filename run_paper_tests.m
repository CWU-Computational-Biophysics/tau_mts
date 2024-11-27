% run_paper_tests

% clear the workspace
clear;

% define save directory
save_dir = fullfile('data', 'nov26_paper_test_run_true');
[~, ~, ~] = mkdir(save_dir);

% calculate map6 values
tau_map6_ratio = [10];
default_params = genParams();
map6_vals = default_params.tau_on .* (1./tau_map6_ratio);

% iterate over values of map_on
for map6_val = map6_vals
    % repeat the simulation five times
    for i = 1:5
        % generate parameters
        params = genParams( ...
            ttot=60*10, ...
            tau_gating=true, ...
            map6_on=map6_val*10, ...
            map6_off=map6_val*10, ...
            tau_on=default_params.tau_on*10, ...
            tau_off=default_params.tau_off*10 ...
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
