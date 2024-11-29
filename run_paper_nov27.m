% run_paper_nov27.m

% producing four figures for the final paper data
% see rough-params.xlsx for more param details

% clear workspace
clear;

% define the number of runs for each parameter
runs = 1:5;

% define static params
dt = 0.01;
ttot = 1000;

% define data directory
save_dir = fullfile('paper_data2');
[~, ~, ~] = mkdir(save_dir);

% define the list of vars to clear after saving data
clear_vars_list = {'mt_grid', 'mt_length', 'mt_state', 'mt_grids', 'growths'};

%{
% fig 1 data
% tau_on = vary, tau_off = 25
% map6_on = 0.25, map6_off = 0.25

% define fig1 save directory
fig1_dir = fullfile(save_dir, 'fig1');
[~, ~, ~] = mkdir(fig1_dir);
fprintf("\nStarting simulations for fig1 @ '%s'\n", fig1_dir);

% true
true_dir = fullfile(fig1_dir, 'fig1_true');
[~, ~, ~] = mkdir(true_dir);
% select a value for tau_on
for tau_val = [25, 50, 100, 200, 250] 
    % repeat 'runs' times
    for i = runs
        % generate parameters
        params = genParams( ...
            ttot=ttot, ...
            dt=dt, ...
            tau_gating=true, ...
            tau_on=tau_val, ...
            tau_off=25, ...
            map6_on=0.25, ...
            map6_off=0.25);

        % print banner
        fprintf("\nSim %i @ tau_on=%g, tau_gating=%s\n", i, params.tau_on, mat2str(params.tau_gating));

        % run the simulation
        try
            [mt_grid, mt_length, mt_state, grids, growths] = simTauMap6(params);

            % save the simulation
            save_name = sprintf('true_sim_%g_%i.mat', params.tau_on, i);
            save_file = fullfile(true_dir, save_name);
            % save(save_file, '-nocompression', '-v7');
            save(save_file, '-v7.3');
            fprintf("Simulation saved to '%s'\n", save_file);

            % delete large vars
            clear(clear_vars_list{:})
        catch exception
            warning("Exception, skipping\n")
        end
    end
end

% false
false_dir = fullfile(fig1_dir, 'fig1_false');
[~, ~, ~] = mkdir(false_dir);
% select a value for tau_on
for tau_val = [25, 50, 100, 200, 250]
    % repeat 'runs' times
    for i = runs
        % generate parameters
        params = genParams( ...
            ttot=ttot, ...
            dt=dt, ...
            tau_gating=false, ...
            tau_on=tau_val, ...
            tau_off=25, ...
            map6_on=0.25, ...
            map6_off=0.25);

        % print banner
        fprintf("\nSim %i @ tau_on=%g, tau_gating=%s\n", i, params.tau_on, mat2str(params.tau_gating));

        % run the simulation
        try
            [mt_grid, mt_length, mt_state, grids, growths] = simTauMap6(params);
    
            % save the simulation
            save_name = sprintf('false_sim_%g_%i.mat', params.tau_on, i);
            save_file = fullfile(false_dir, save_name);
            % save(save_file, '-nocompression', '-v7');
            save(save_file, '-v7.3');
            fprintf("Simulation saved to '%s'\n", save_file);

            % delete large vars
            clear(clear_vars_list{:})
        catch exception
            warning("Exception, skipping\n")
        end
    end
end
%}

% fig 2 data
% tau_on = vary, tau_off = 0.25
% map6_on = 0.25, map6_off = 0.25

% define fig2 save directory
fig2_dir = fullfile(save_dir, 'fig2');
[~, ~, ~] = mkdir(fig2_dir);
fprintf("\nStarting simulations for fig2 @ '%s'\n", fig2_dir);

% true
true_dir = fullfile(fig2_dir, 'fig2_true');
[~, ~, ~] = mkdir(true_dir);
% select a value for tau_on
for tau_val = [0.25, 1, 2.5, 5, 10, 20, 25]
    % repeat 'runs' times
    for i = runs
        % generate parameters
        params = genParams( ...
            ttot=ttot, ...
            dt=dt, ...
            tau_gating=true, ...
            tau_on=tau_val, ...
            tau_off=0.25, ...
            map6_on=0.25, ...
            map6_off=0.25);

        % print banner
        fprintf("\nSim %i @ tau_on=%g, tau_gating=%s\n", i, params.tau_on, mat2str(params.tau_gating));

        % run the simulation
        try
            [mt_grid, mt_length, mt_state, grids, growths] = simTauMap6(params);
    
            % save the simulation
            save_name = sprintf('true_sim_%g_%i.mat', params.tau_on, i);
            save_file = fullfile(true_dir, save_name);
            % save(save_file, '-nocompression', '-v7');
            save(save_file, '-v7.3');
            fprintf("Simulation saved to '%s'\n", save_file);

            % delete large vars
            clear(clear_vars_list{:})
        catch exception
            warning("Exception, skipping\n")
        end
    end
end

% false
false_dir = fullfile(fig2_dir, 'fig2_false');
[~, ~, ~] = mkdir(false_dir);
% select a value for tau_on
for tau_val = [0.25, 1, 2.5, 5, 10, 20, 25]
    % repeat 'runs' times
    for i = runs
        % generate parameters
        params = genParams( ...
            ttot=ttot, ...
            dt=dt, ...
            tau_gating=false, ...
            tau_on=tau_val, ...
            tau_off=0.25, ...
            map6_on=0.25, ...
            map6_off=0.25);

        % print banner
        fprintf("\nSim %i @ tau_on=%g, tau_gating=%s\n", i, params.tau_on, mat2str(params.tau_gating));

        % run the simulation
        try
            [mt_grid, mt_length, mt_state, grids, growths] = simTauMap6(params);
    
            % save the simulation
            save_name = sprintf('false_sim_%g_%i.mat', params.tau_on, i);
            save_file = fullfile(false_dir, save_name);
            % save(save_file, '-nocompression', '-v7');
            save(save_file, '-v7.3');
            fprintf("Simulation saved to '%s'\n", save_file);

            % delete large vars
            clear(clear_vars_list{:})
        catch exception
            warning("Exception, skipping\n")
        end
    end
end

%{
% fig 3 data
% tau_on = 0.25, tau_off = 0.25
% map6_on = vary, map6_off = 0.25

% define fig3 save directory
fig3_dir = fullfile(save_dir, 'fig3');
[~, ~, ~] = mkdir(fig3_dir);
fprintf("\nStarting simulations for fig3 @ '%s'\n", fig3_dir);

% true
true_dir = fullfile(fig3_dir, 'fig3_true');
[~, ~, ~] = mkdir(true_dir);
% select a value for tau_on
for map6_val = [0.25, 1, 2.5, 5, 10, 20, 25]
    % repeat 'runs' times
    for i = runs
        % generate parameters
        params = genParams( ...
            ttot=ttot, ...
            dt=dt, ...
            tau_gating=true, ...
            tau_on=0.25, ...
            tau_off=0.25, ...
            map6_on=map6_val, ...
            map6_off=0.25);

        % print banner
        fprintf("\nSim %i @ map6_on=%g, tau_gating=%s\n", i, params.map6_on, mat2str(params.tau_gating));

        % run the simulation
        try
            [mt_grid, mt_length, mt_state, grids, growths] = simTauMap6(params);
    
            % save the simulation
            save_name = sprintf('true_sim_%g_%i.mat', params.map6_on, i);
            save_file = fullfile(true_dir, save_name);
            % save(save_file, '-nocompression', '-v7');
            save(save_file, '-v7.3');
            fprintf("Simulation saved to '%s'\n", save_file);

            % delete large vars
            clear(clear_vars_list{:})
        catch exception
            warning("Exception, skipping\n")
        end
    end
end

% false
false_dir = fullfile(fig3_dir, 'fig3_false');
[~, ~, ~] = mkdir(false_dir);
% select a value for tau_on
for map6_val = [0.25, 1, 2.5, 5, 10, 20, 25]
    % repeat 'runs' times
    for i = runs
        % generate parameters
        params = genParams( ...
            ttot=ttot, ...
            dt=dt, ...
            tau_gating=false, ...
            tau_on=0.25, ...
            tau_off=0.25, ...
            map6_on=map6_val, ...
            map6_off=0.25);

        % print banner
        fprintf("\nSim %i @ map6_on=%g, tau_gating=%s\n", i, params.map6_on, mat2str(params.tau_gating));

        % run the simulation
        try
            [mt_grid, mt_length, mt_state, grids, growths] = simTauMap6(params);
    
            % save the simulation
            save_name = sprintf('false_sim_%g_%i.mat', params.map6_on, i);
            save_file = fullfile(false_dir, save_name);
            % save(save_file, '-nocompression', '-v7');
            save(save_file, '-v7.3');
            fprintf("Simulation saved to '%s'\n", save_file);

            % delete large vars
            clear(clear_vars_list{:})
        catch exception
            warning("Exception, skipping\n")
        end
    end
end


% fig 4 data
% tau_on = 0.25, tau_off = 0.25
% map6_on = vary, map6_off = 0.25

% define fig4 save directory
fig4_dir = fullfile(save_dir, 'fig4');
[~, ~, ~] = mkdir(fig4_dir);
fprintf("\nStarting simulations for fig4 @ '%s'\n", fig4_dir);

% true
true_dir = fullfile(fig4_dir, 'fig4_true');
[~, ~, ~] = mkdir(true_dir);
% select a value for tau_on
for map6_val = [0.25, 1, 2.5, 5, 10, 20, 25]
    % repeat 'runs' times
    for i = runs
        % generate parameters
        params = genParams( ...
            ttot=ttot, ...
            dt=dt, ...
            tau_gating=true, ...
            tau_on=0.25, ...
            tau_off=0.25, ...
            map6_on=map6_val, ...
            map6_off=0.25);

        % print banner
        fprintf("\nSim %i @ map6_on=%g, tau_gating=%s\n", i, params.map6_on, mat2str(params.tau_gating));

        % run the simulation
        try
            [mt_grid, mt_length, mt_state, grids, growths] = simTauMap6(params);
    
            % save the simulation
            save_name = sprintf('true_sim_%g_%i.mat', params.map6_on, i);
            save_file = fullfile(true_dir, save_name);
            % save(save_file, '-nocompression', '-v7');
            save(save_file, '-v7.3');
            fprintf("Simulation saved to '%s'\n", save_file);

            % delete large vars
            clear(clear_vars_list{:})
        catch exception
            warning("Exception, skipping\n")
        end
    end
end

% false
false_dir = fullfile(fig4_dir, 'fig4_false');
[~, ~, ~] = mkdir(false_dir);
% select a value for tau_on
for map6_val = [0.25, 1, 2.5, 5, 10, 20, 25]
    % repeat 'runs' times
    for i = runs
        % generate parameters
        params = genParams( ...
            ttot=ttot, ...
            dt=dt, ...
            tau_gating=false, ...
            tau_on=0.25, ...
            tau_off=0.25, ...
            map6_on=map6_val, ...
            map6_off=0.25);

        % print banner
        fprintf("\nSim %i @ map6_on=%g, tau_gating=%s\n", i, params.map6_on, mat2str(params.tau_gating));

        % run the simulation
        try
            [mt_grid, mt_length, mt_state, grids, growths] = simTauMap6(params);
    
            % save the simulation
            save_name = sprintf('false_sim_%g_%i.mat', params.map6_on, i);
            save_file = fullfile(false_dir, save_name);
            % save(save_file, '-nocompression', '-v7');
            save(save_file, '-v7.3');
            fprintf("Simulation saved to '%s'\n", save_file);

            % delete large vars
            clear(clear_vars_list{:})
        catch exception
            warning("Exception, skipping\n")
        end
    end
end
%}
