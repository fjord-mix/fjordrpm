function [p, a, f, t] = example_intermediate_circulation

% EXAMPLE_INTERMEDIATE_CIRCULATION  Example to run the box model with intermediate circulation.
%   P = EXAMPLE_INTERMEDIATE_CIRCULATION sets the user-defined box model parameters
%   P for the intermediate circulation example.

[p,~] = get_model_default_parameters();

%% Fjord geometry parameters
p.N = 6; % number of above-sill model layers
p.sill = 0; % flag for sill (0 = no sill, 1 = sill)

%% Tuning/artificial parameters
p.Hmin = NaN; % minimum box thickness, NaN = no min thickness
p.Snudge = 0.*ones(1, p.N-1); % controls layer thickness for nudging, 0 if no nudging

%% Time step
p.dt = 0.1; % time stepping (units are days)
p.t_end = 10; % time to end the simulation
t = 0:p.dt:p.t_end;

%% Forcing parameters: idealised forcing
% Shelf forcing parameters
p.z0 = 50; % sets strength of shelf stratification
p.zd = 30; % sets strength of shelf oscillation, zero if no oscillation in shelf
p.tw = 10; % oscillation period of shelf forcing (days)
p.sf = @(S1, S2, Z, Z0) S1 - (S1 - S2)*exp(Z./Z0); % functional form of the shelf stratification
p.Stop = 30; % shelf salinity at the ocean surface
p.Sbottom = 35; % shelf salinity at the ocean floor
p.Ttop = 3; % shelf temperature at the ocean surface
p.Tbottom = 3; % shelf temperature at the ocean floor
% Glacier forcing parameters
p.Qv0 = 0; % volumetric flow rate of subglacial discharge
p.P0 = 0; % plume entrainment width, 0 = no plume
p.gf = @(T, Q) Q*ones(size(T)); % functional form of subglacial discharge rate depending on time
% Iceberg forcing parameters
p.M0 = 0; % iceberg melt efficiency, 0 = no icebergs
p.nu0 = 25; % iceberg volume profile coefficient
p.E0 = 1e-7; % iceberg export efficiency
p.if = @(NU, H, Z) (NU/H)*exp(NU*Z/H)/(1-exp(-NU)); % functional form of iceberg depth profile    

% forcing structure
f = get_idealised_forcing(p, t);

%% Initial conditions
a = get_initial_conditions(p, f);

end