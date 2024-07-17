% tau_map6_run.m
% a standard script for sweeping parameters

% iterate over parameter sets
% for i = [1, 5, 10, 50]
% 	for j = [true, false]
% 		for k = [1]
% 			% load params
% 			tau_map6_params;
% 
% 			% update values
% 			tm_ratio = i;
% 			t_force = j;
% 			ttot = 100;
% 
% 			% update save
% 			data_dir = sprintf("data/paper_anims" + ttot + "_" + t_force);
% 			export_dir = fullfile(data_dir);
% 			sim_name = sprintf("taumap6_" + tm_ratio + "_" + t_force + "_" + k);
% 
% 			% run simulation
% 			tau_map6;
% 		end
% 	end
% end


% load params
tau_map6_params;

% update values
tm_ratio = 50;
t_force = true;
ttot = 500;

% update save
data_dir = sprintf("data/paper_anims_test" + ttot + "_" + t_force);
export_dir = fullfile(data_dir);
sim_name = sprintf("taumap6_" + tm_ratio + "_" + t_force + "_" + k);

% run simulation
tau_map6;