% EXAMPLE_ICEBERGS  Run the zmodel with example iceberg parameters.
addpath(genpath('./'))
clearvars
close all
% Get default parameters which can be overwritten for this specific
% example.
[p, t, f, a] = get_model_default_parameters;

%% Set the input parameters p.
% Set ZMODEL layer properties.
p.N = 80; % number of above-sill model layers

% Set fjord geometry parameters.
p.sill = 1; % flag for sill (0 = no sill, 1 = sill)

% Set shelf oscillation to zero.
p.zd = 0; % sets strength of shelf oscillation, zero if no oscillation in shelf

% Set plume forcing to zero.
p.Qv0 = 0; % volumetric flow rate of subglacial discharge

% Set iceberg forcing.
p.icestatic = 1; % 1 for static icebergs, 0 for evolving
p.M0 = 7e-7; % iceberg melt efficiency
p.A0 = 2e9; % scaling of iceberg area
p.U0 = 0.25; % scale upwelling
p.nu0 = 25; % iceberg profile coefficient
% p.E0 = 1e-7; % iceberg export efficiency p.uIce = 0.005; % iceberg
% down-fjord velocity
p.if = @(NU, H, Z) (NU/H)*exp(NU*Z/H)/(1-exp(-NU)); % functional form of iceberg depth profile

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

name = 'example_icebergs';
% Save the fjord structure, including input parameters, initial conditions,
% boundary conditions and solution.
save([output_folder, name, '.mat'], 'p', 't', 'f', 'a', 's');
