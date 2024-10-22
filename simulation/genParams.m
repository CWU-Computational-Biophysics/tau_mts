function params = genParams(varargin)
%GENPARAMS Summary of this function goes here
%   Detailed explanation goes here

% initialize the input parser
input = inputParser;

% define input validation functions
validScalar = @(x) isscalar(x) && (x >= 0);
validInt = @(x) validScalar(x) && (mod(x, 1) == 0);

% simulation input parameters
% [second]
addParameter(input, 'dt', 0.01, validScalar);
addParameter(input, 'ttot', 1, validScalar);

% constrained input parameters
% [#]
addParameter(input, 'init_domains', 20, validInt);
addParameter(input, 'alpha_t', 0.0, validScalar);
addParameter(input, 'alpha_m', 0.0, validScalar);
% [micrometer]
addParameter(input, 'dx', 0.036, validScalar);
% [1 / second]
addParameter(input, 'f_res', 50, validScalar);
addParameter(input, 'f_cat', 30, validScalar);
% [mirometer / second]
addParameter(input, 'vp', 3.6, validScalar);
addParameter(input, 'vm', 3.6, validScalar);

% tunable input parameters
% [logical]
addParameter(input, 'tau_gating', true, @islogical);
% [1 / second]
addParameter(input, 'map6_on', 0.1, validScalar);
addParameter(input, 'map6_off', 0.1, validScalar);
addParameter(input, 'tau_on', 0.1, validScalar);
addParameter(input, 'tau_off', 0.1, validScalar);

% parse inputs and output to struct
parse(input, varargin{:});
params = input.Results;

end