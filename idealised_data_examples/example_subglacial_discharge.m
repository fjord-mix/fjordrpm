% EXAMPLE_SUBGLACIAL_DISCHARGE  Run the zmodel with example
% subglacial discharge parameters.
addpath(genpath('./'))
clearvars
close all
% Get default parameters which can be overwritten for this specific
% example.
[p, t, f, a] = get_model_default_parameters;

%% Set the input parameters p.
% ZMODEL layer properties
p.N = 60; % number of above-sill model layers

% Fjord geometry parameters
p.sill = 0; % flag for sill (0 = no sill, 1 = sill)

% Set shelf oscillation to zero.
p.zd = 0; % sets strength of shelf oscillation, zero if no oscillation in shelf

% Set plume forcing.
p.Qv0 = 100; % volumetric flow rate of subglacial discharge
p.P0 = 25; % plume entrainment width, 0 = no plume
p.gf = @(T, Q) Q*ones(size(T)); % functional form of subglacial discharge rate depending on time

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

name = 'example_subglacial_discharge';
% Save the fjord structure, including input parameters, initial conditions,
% boundary conditions and solution.
save([output_folder, name, '.mat'], 'p', 't', 'f', 'a', 's');
