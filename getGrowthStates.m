function vals = getGrowthStates()
%GETGROWTHSTATES Summary of this function goes here
%   Detailed explanation goes here

% define the growth state struct
vals = struct();
vals.STABLE = 0;
vals.GROWING = 1;
vals.SHRINKING = -1;

end