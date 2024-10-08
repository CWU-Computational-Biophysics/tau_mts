% tau_map6_run.m
% a standard script for sweeping parameters

% clear workspace
clear;

% iterate over parameter sets
for i = [0.1, 1, 5, 10, 50]
	for j = [true, false]
		for t = [50, 100]
			for k = [1]
				% load params
				tau_map6_params;

				% update values
				tm_ratio = i;
				t_force = j;
				ttot = t;

				% update save
				data_dir = sprintf("data" + "/mapisdot1/paper_" + ttot + "_" + t_force);
				export_dir = fullfile(data_dir);
				sim_name = sprintf("taumap6_" + round(tm_ratio) + "_" + t_force + "_" + k);

				% run simulation
				tau_map6;
			end
		end
	end
end
