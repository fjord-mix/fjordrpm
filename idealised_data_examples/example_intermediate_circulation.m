% EXAMPLE_INTERMEDIATE_CIRCULATION  Run the zmodel with example
% intermediate circulation parameters.
addpath(genpath('./'))
clearvars
close all
% Get default parameters which can be overwritten for this specific
% example.
[p, t, f, a] = get_model_default_parameters;

%% Set the input parameters p.
% ZMODEL layer properties.
p.N = 60; % number of above-sill model layers

% Fjord geometry parameters
p.sill = 0; % flag for sill (0 = no sill, 1 = sill)

% Set shelf forcing. 
p.z0 = 50; % sets strength of shelf stratification
p.zd = 30; % sets strength of shelf oscillation, zero if no oscillation in shelf
p.tw = 10; % oscillation period of shelf forcing (days)
p.sf = @(S1, S2, Z, Z0) S1 - (S1 - S2)*exp(Z./Z0); % functional form of the shelf stratification
p.Stop = 30; % shelf salinity at the ocean surface
p.Sbottom = 35; % shelf salinity at the ocean floor
p.Ttop = 3; % shelf temperature at the ocean surface
p.Tbottom = 3; % shelf temperature at the ocean floor

% Set plume forcing to zero.
p.Qv0 = 0; % volumetric flow rate of subglacial discharge

% Set iceberg forcing to zero.
p.M0 = 0; % iceberg melt efficiency, 0 = no icebergs

%% Set the input time vector t.
% Time values at which to compute the solution (in days)
t = 0:1:200;

%% Set the boundary conditions f.
% Boundary conditions for the input parameters at each timestep
f = get_idealised_forcing(p, t);

%% Set the initial conditions a. 
% Initial conditions for the input parameters with the initial boundary
% conditions.
a = get_initial_conditions(p, f);

%% Run the ZMODEL. 
s = zmodel(p, t, f, a);

%% Save the results.
% Choose where to save the results and make the directory if it doesn't
% exist.
output_folder='./idealised_data_examples/results/';
if not(isfolder(output_folder))
    mkdir(output_folder)
end

name = 'example_intermediate_circulation';
% Save the fjord structure, including input parameters, initial conditions,
% boundary conditions and solution.
save([output_folder, name, '.mat'], 'p', 't', 'f', 'a', 's');
