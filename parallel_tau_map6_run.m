% tau_map6_run.m
% a standard script for sweeping parameters

% clear workspace
clear;

% Define parameter sets
i_values = [1, 5, 10, 50];
j_values = [true, false];
t_values = [50, 100];
k_values = [1, 2, 3];

% Create a parallel pool
parpool;

% Iterate over parameter sets using parfor
parfor idx = 1:length(i_values) * length(j_values) * length(t_values) * length(k_values)
    % Calculate indices for each parameter
    [i_idx, j_idx, t_idx, k_idx] = ind2sub([length(i_values), length(j_values), length(t_values), length(k_values)], idx);

    % Load params
    tau_map6_params;

    % Update values
    tm_ratio = i_values(i_idx);
    t_force = j_values(j_idx);
    ttot = t_values(t_idx);
    k = k_values(k_idx);

    % Update save
    data_dir = sprintf("data/paper_%d_%d", ttot, t_force);
    export_dir = fullfile(data_dir);
    sim_name = sprintf("taumap6_%d_%d_%d", tm_ratio, t_force, k);

    % Run simulation
    tau_map6;
end

% Shut down the parallel pool
delete(gcp('nocreate'));
