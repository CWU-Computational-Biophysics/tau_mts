% tau_map6.m
% Tau/Map6 project coarse grain model

% load user-input parameters
% clear;
% tau_map6_params;
% tau_map6_params_default;

% calculated parameters
% base binding rate for tau [1/second] (T0)
t_on = tm_ratio * m_on;
% unbinding rate for tau [1/second] (Toff)
t_off = t_on;
% total number of time steps [#] (N)
time_scale = 10;
ttot_prime = ttot/time_scale;
steps = (ttot_prime/dt) + 1;
time = (dt*time_scale) .* (1:steps);

% initialize the MT grid array (MTgrid)
mt_grid = ones(steps, grid_points_init);
% relate the grid state to an integer
state_notexist = -1;
state_empty = 0;
state_tau = 1;
state_map6 = 2;
% empty the grid
mt_grid = mt_grid .* state_notexist;
% empty the first row
mt_grid(1, 1:grid_points_init) = state_empty;

% initialize the MT length array (MTlength)
% 1: length [#],
% 2: growthstate [#],
% 3: taufractip [#],
% 4: taufraclength [#],
% 5: tauplusendasym [#],
% 6: mapfractip [#],
% 7: mapplusendasym [#]
mt_length = zeros(steps, 2);
mt_length(1,1) = grid_points_init;
% relate the growth state to an integer
state_static = 0;
state_shrinking = -1;
state_growing = 1;
mt_length(1,2) = state_static;

% parameter readout
fprintf("sim name: '" + sim_name + "'\n");
fprintf("tm_ratio: " + tm_ratio + "\n");
fprintf("t_force: " + t_force + "\n");
fprintf("sim time: " + ttot + " seconds\n");


% simulation begin
% for loop over time
for step_i = 1:(steps-1)
	% get the current number of grid points
	[~, grid_points] = find(mt_grid(step_i, :)==state_notexist, 1);

	% check if grid_points array is empty (mt fully exists)
	if numel(grid_points) == 0
		grid_points = numel(mt_grid(step_i, :));
	else
		grid_points = grid_points-1;
	end

	% pass on existing sites
	mt_grid(step_i+1, 1:grid_points) = state_empty;

	% loop over grid points
	for grid_i = 1:grid_points
		% determine tau onrate based on next neighbor occupancy
		tau_next = 0;
		map6_next = 0;

		% check left protein
		if grid_i ~= 1
			left_protein = mt_grid(step_i, grid_i-1);
			if left_protein == state_tau
				tau_next = tau_next + 1;
			elseif left_protein == state_map6
				map6_next = map6_next + 1;
			end
		end

		% check right protein
		if grid_i ~= grid_points
			right_protein = mt_grid(step_i, grid_i+1);
			if right_protein == state_tau
				tau_next = tau_next + 1;
			elseif right_protein == state_map6
				map6_next = map6_next + 1;
			end
		end

		% calculate effective offrate
		t_off_eff = t_off - alpha_t*t_off*tau_next;
		m_off_eff = m_off - alpha_t*m_off*map6_next;

		% update occupancy numbers
		site_protein = mt_grid(step_i, grid_i);

		% empty site
		if site_protein == state_empty
			if rand() < t_on*dt
				% tau binding
				mt_grid(step_i+1, grid_i) = state_tau;
			elseif rand() < m_on*dt
				% map6 binding
				mt_grid(step_i+1, grid_i) = state_map6;
			else
				% no change
				mt_grid(step_i+1, grid_i) = site_protein;
			end
		end

		% tau site
		if site_protein == state_tau
			if rand() < t_off_eff*dt
				% tau unbinding
				mt_grid(step_i+1, grid_i) = state_empty;
			else
				% no change
				mt_grid(step_i+1, grid_i) = site_protein;
			end
		end

		% map6 site
		if site_protein == state_map6
			if rand() < m_off_eff*dt
				% map6 unbinding
				mt_grid(step_i+1, grid_i) = state_empty;
			else
				% no change
				mt_grid(step_i+1, grid_i) = site_protein;
			end
		end
	end

	% chance for MT state change
	growth_state = mt_length(step_i, 2);

	% dynamic MTs when tau bound in last three sites
	tau_sum = 0;
	if grid_points >= 3
		tau_grid = mt_grid(step_i, (grid_points-2):grid_points)==state_tau;
		tau_sum = sum(tau_grid);
	end

	if tau_sum > 0
		% chance to change growth state
		if growth_state == state_static
			growth_state = state_growing;
		elseif growth_state == state_shrinking
			if rand() < f_res*dt
				growth_state = state_growing;
			end
		elseif growth_state == state_growing
			if rand() < f_cat*dt
				growth_state = state_shrinking;
			end
		end
	else
		% MT is not dynamic
		growth_state = 0;
	end

	% enact growing or shrinking
	if growth_state == state_growing
		% increase the number of rows
		% first define a new column of the not exist state
		new_col = ones(steps, 1);
		new_col = new_col .* state_empty;
		
		% append the column to the array
		mt_grid = [mt_grid, new_col];
		
		% update the number of grid points
		grid_points = grid_points + 1;

		% empty the grid at the next step
		mt_grid(step_i+1, grid_points) = state_empty;

		% check for force tau
		if t_force
			mt_grid(step_i+1, grid_points) = state_tau;
		end
	elseif growth_state == state_shrinking
		% update the number of grid points
		grid_points = grid_points - 1;
	end

	% get column count
	col_count = numel(mt_grid(step_i+1, :));

	% check if any grid points need to be marked non existent
	if col_count > grid_points
		% mark non existent between grid points (not-inclusive) and max entry
		mt_grid(step_i+1, (grid_points+1):end) = state_notexist;
	end

	% pass on grid point count and growth state
	mt_length(step_i+1, 1) = grid_points;
	mt_length(step_i+1, 2) = growth_state;

	% check for MT death
	if grid_points == 0
		break
	end
end


% get array states
final_n = 10;
if grid_points > final_n
	tip_array = mt_grid(end, (grid_points-final_n):grid_points);
	length_array = mt_grid(end, 1:(grid_points-(final_n+1)));

	% calculate final binding fraction asymmetry
	tau_frac_tip = sum(tip_array==state_tau) / numel(tip_array);
	tau_frac_length = sum(length_array==state_tau) / numel(length_array);
	tau_plus_end_asym = tau_frac_tip / tau_frac_length;
	map6_frac_tip = sum(tip_array==state_map6) / numel(length_array);
	map6_frac_length = sum(length_array==state_map6) / numel(length_array);
	map6_plus_end_asym = map6_frac_tip / map6_frac_length;
end

% export the workspace in python compatible .mat
export_file = fullfile(export_dir, sim_name);
[~, ~, ~] = mkdir(export_dir);
save(export_file, "-nocompression", "-v7");
fprintf("sim saved @ '" + export_file + "'.\n");
