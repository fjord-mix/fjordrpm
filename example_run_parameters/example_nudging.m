function run_output = example_nudging

% EXAMPLE_NUDGING  Example to run the box model with nudging.
%   P = EXAMPLE_NUDGING sets the user-defined box model parameters
%   P for the nudging example.

[p,~] = get_model_default_parameters();

%% Fjord geometry parameters
p.N         = 6; % number of above-sill model layers
p.sill      = 0; % flag for sill (0 = no sill, 1 = sill)

%% Tuning/artificial parameters
p.Hmin   = 0; % minimum box thickness, NaN = no min thickness
p.trelax = 100; % controls layer nudging, NaN = no nudging
p.Snudge = [32,33,34,34.5,34.8].*ones(1, p.N-1); % controls layer thickness for nudging, 0 if no nudging

%% Time step
p.dt    = 0.1; % time stepping (units are days)
p.t_end = 100; % time to end the simulation
t       = 0:p.dt:p.t_end;

%% Forcing parameters: idealised forcing
% Shelf forcing parameters
p.z0 = 50; % sets strength of shelf stratification
p.zd = 0; % sets strength of shelf oscillation, zero if no oscillation in shelf
p.tw = 0; % oscillation period of shelf forcing (days)
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
run_output.f = get_idealised_forcing(p, t);

%% Initial conditions
run_output.a = get_initial_conditions(p, run_output.f);

%% adding the other structures to output structure
run_output.p = p;
run_output.t = t;

end