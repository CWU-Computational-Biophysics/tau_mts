function [mt_grid, mt_length, mt_state, grids, growths] = simTauMap6(params)
%SIMTAUMAP6 Summary of this function goes here
%   Detailed explanation goes here

%% Load and validate inputs
% load params via input parser
input = inputParser;
addRequired(input, 'params', @(x) isstruct(x) && ~isempty(x));
parse(input, params);
params = input.Results.params;

% validate that ttot is an integer multiple of dt
if mod(params.ttot, params.dt) ~= 0
    error('ttot (%g) must be an integer multiple of dt (%g).', params.ttot, params.dt)
end

% validate that rate*dt < 1 and 0.1
rate_vector = [params.f_res, params.f_cat, params.map6_on, params.map6_off, params.tau_on, params.tau_off];
eff_rate_vector = rate_vector .* params.dt;
if any(eff_rate_vector >= 1)
    error('dt=%g is too high, some rate parameter is >= 1', params.dt)
elseif any(eff_rate_vector >= 0.1)
    warning('dt=%g is high, some rate parameter is >= 0.1', params.dt)
end

% validate that effective growth rate is below 1 domain/step
eff_growth_rate = floor(params.vp * params.dt / params.dx);
if eff_growth_rate > 1
    warning('choice of dt (%g) is high, growth rate is multiple domains per step (%g)', params.dt, eff_growth_rate)
end

%% Initialize the simulation
% get the simulation fields
grids = getGridStates();
growths = getGrowthStates();

% calculate time steps [#]
nsteps = (params.ttot / params.dt) + 1;
% time = params.dt .* (1:nsteps);

% initialize the MT grid array with nonexistent grid states
mt_grid = ones(nsteps, params.init_domains) .* grids.NOTEXIST;
% empty the first row
mt_grid(1, 1:params.init_domains) = grids.EMPTY;

% initialize the MT length array
mt_length = zeros(nsteps, 1);
mt_length(1) = params.init_domains * params.dx;

% initialize the MT state array
mt_state = zeros(nsteps, 1);
mt_state(1) = growths.STABLE;


%% Simulation start
for si = 1:(nsteps-1)
    %% Current grid vector stuff
    % calculate the current number of grid points from the MT length
    % this method locates the first index of a NOTEXIST grid point
    curr_grid_vector = mt_grid(si, :);
    [~, curr_gps] = find(curr_grid_vector==grids.NOTEXIST, 1);

    % check if curr_gps is empty (MT fully exists / no NOTEXIST spaces)
    % else, the number of grid spaces in existence is one less
    if numel(curr_gps) == 0
        curr_gps = numel(mt_grid(si, :));
    else
        curr_gps = curr_gps - 1;
    end

    % pass the current grid vector to the next step as empty
    mt_grid(si+1, 1:curr_gps) = grids.EMPTY;

    %% Individual grid point checks
    % loop ovver grid points
    for gi = 1:curr_gps
        % calculate neighbor occupancy
        tau_next = 0;
        map6_next = 0;

        % check left protein
        if gi ~= 1
            left_protein = mt_grid(si, gi-1);
            if left_protein == grids.TAU
                tau_next = tau_next + 1;
            elseif left_protein == grids.MAP6
                map6_next = map6_next + 1;
            end
        end

        % check right protein
        if gi ~= curr_gps
            right_protein = mt_grid(si, gi+1);
            if right_protein == grids.TAU
                tau_next = tau_next + 1;
            elseif right_protein == grids.MAP6
                map6_next = map6_next + 1;
            end
        end

        % calculate protein effective on/off rates
        tau_off_eff = (params.tau_off - (params.alpha_t * params.tau_off * tau_next)) * params.dt;
        map6_off_eff = (params.map6_off - (params.alpha_m * params.map6_off * map6_next)) * params.dt;
        tau_on_eff = params.tau_on * params.dt;
        map6_on_eff = params.map6_on * params.dt;

        % get the current protein for occupancy update
        curr_protein = mt_grid(si, gi);

        % empty site options (new protein binding)
        if curr_protein == grids.EMPTY
            if rand() < tau_on_eff
                % tau binding
                mt_grid(si+1, gi) = grids.TAU;
            elseif rand() < map6_on_eff
                % map6 binding
                mt_grid(si+1, gi) = grids.MAP6;
            else
                % no change
                mt_grid(si+1, gi) = grids.EMPTY;
            end
        end

        % tau site options (unbinding or not)
        if curr_protein == grids.TAU
            if rand() < tau_off_eff
                % tau unbinding
                mt_grid(si+1, gi) = grids.EMPTY;
            else
                % no change
                mt_grid(si+1, gi) = grids.TAU;
            end
        end

        % map6 site options (unbinding or not)
        if curr_protein == grids.MAP6
            if rand() < map6_off_eff
                % map6 unbinding
                mt_grid(si+1, gi) = grids.EMPTY;
            else
                % no change
                mt_grid(si+1, gi) = grids.MAP6;
            end
        end
    end


    %% Check for MT dynamics
    % define the boolean grid for tau proteins
    % this grid uses the last three existing sites
    last = curr_gps;
    first = max([curr_gps - 2, 1]);
    tau_grid = mt_grid(si, (first:last))==grids.TAU;
    tau_sum = sum(tau_grid);

    % MT may become dynamic when tau bound in the last three sites
    curr_state = mt_state(si);
    next_state = mt_state(si);
    if tau_sum > 0
        % chance to change to growth state
        if curr_state == growths.STABLE
            % always go from stable to growing
            next_state = growths.GROWING;
        elseif curr_state == growths.SHRINKING
            % chance for rescue
            if rand() < (params.f_res * params.dt)
                next_state = growths.GROWING;
            end
        elseif curr_state == growths.GROWING
            % chance for catastrophe
            if rand() < (params.f_cat * params.dt)
                next_state = growths.SHRINKING;
            end
        end
    else
        next_state = growths.STABLE;
    end
    mt_state(si+1) = next_state;

    %% Check for MT growth/shrinking
    % update MT length
    if curr_state == growths.GROWING
        % increase the MT length
        mt_length(si+1) = mt_length(si) + (params.vp * params.dt);
    elseif curr_state == growths.SHRINKING
        % decrease the MT length
        mt_length(si+1) = mt_length(si) - (params.vm * params.dt);
    else
        % persist MT length
        mt_length(si+1) = mt_length(si);
    end

    %% Update grid points based on growth
    % calculate the next number of grid points from the MT length
    length_based_grid_count = floor(mt_length(si+1) / params.dx);

    % calculate the next number of grid points from the mt grid
    % this is done by finding the first index of NOTEXIST
    next_grid_vector = mt_grid(si+1, :);
    next_step_grid_max = length(next_grid_vector);
    [~, next_step_grid_count] = find(next_grid_vector==grids.NOTEXIST, 1);

    % check if next_gps is empty (MT fully exists / no NOTEXIST spaces)
    % else, the number of grid spaces in existence is one less
    if numel(next_step_grid_count) == 0
        next_step_grid_count = length(next_grid_vector);
    else
        next_step_grid_count = next_step_grid_count - 1;
    end


    % find the difference in grid points
    grid_count_difference = length_based_grid_count - next_step_grid_count;

    % update the mt_grid based on diff_gps
    if grid_count_difference > 0
        % calculate how many grid columns need to be appended (if any)
        appending_grid_count = length_based_grid_count - next_step_grid_max;

        % adding new grid spaces to mt_grid
        if appending_grid_count > 0
            % define a new NOTEXIST column to append to the array
            new_col = ones(nsteps, grid_count_difference) .* grids.NOTEXIST;

            % append the NOTEXIST column
            mt_grid = [mt_grid, new_col]; %#ok<AGROW>
        end

        % empty the new grid points at the next step
        mt_grid(si+1, next_step_grid_count:length_based_grid_count) = grids.EMPTY;

        % check for tau gating
        if params.tau_gating
            mt_grid(si+1, next_step_grid_count:length_based_grid_count) = grids.TAU;
        end
    elseif grid_count_difference < 0
        % removing grid spaces from mt_grid
        % set the diff_gps values to NOTEXIST
        mt_grid(si+1, length_based_grid_count:next_step_grid_count) = grids.NOTEXIST;
    end

    %% Update the grid to reflect potentially nonexistent sites
    % get the column count
    col_count = length(mt_grid(si+1, :));

    % check if any grid points need to be marked NOTEXIST
    if col_count > length_based_grid_count
        % mark the points between true and col as NOTEXIST
        mt_grid(si+1, length_based_grid_count:end) = grids.NOTEXIST;
    end
end
end
