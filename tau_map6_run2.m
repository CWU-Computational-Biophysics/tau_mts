% tau_map6_run2.m
% testing dt invariance

% clear workspace
clear;

% define save location
export_dir = fullfile(sprintf('data/dt_invariance'));

% iterate over three legal values of 'dt'
for dt_val = [0.01, 0.001, 0.0001]
    % load params
    tau_map6_params;

    % update values
    tm_ratio = 1;
    t_force = true;
    ttot = 50;
    dt = dt_val;

    % run each simulation three times
    for i = 1:3
        % update save file name
        sim_name = sprintf('taumap6_%g_%i.mat', dt_val, i);

        % run simulation
        tau_map6;
    end
end