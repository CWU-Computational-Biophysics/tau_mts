% tau_map6_run2.m
% testing dt invariance

% clear workspace
clear;

% iterate over three legal values of 'dt'
for dt_val = [0.01, 0.005, 0.001]
    % load params
    tau_map6_params;

    % define save location
    export_dir = fullfile(sprintf('data/tests3'));

    % update values
    tm_ratio = 1;
    t_force = true;
    ttot = 10;
    dt = dt_val;
    fprintf('dt=%g', dt)

    % run each simulation three times
    for i = 1:3
        % update save file name
        sim_name = sprintf('taumap6_%g_%i.mat', dt_val, i);

        % run simulation
        tau_map6;
    end
end