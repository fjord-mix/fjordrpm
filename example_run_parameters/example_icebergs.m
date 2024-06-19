function run_output = example_icebergs

% EXAMPLE_ICEBERGS  Example to run the zmodel with icebergs.
%   RUN_OUTPUT = EXAMPLE_ICEBERGS returns a RUN_OUTPUT structure containing
%   the user-defined zmodel parameters P and T for the icebergs example,
%   along with the initial conditions A and boundary conditions F.

%% Set the model default parameters.
[p,~] = get_model_default_parameters();

%% Set the specific parameter for this example.
% Set a vector of time values at which to compute the solution (in days).
t       = 0:1:200;

% Set the ZMODEL layer properties.
p.N         = 80; % number of above-sill model layers
p.fixedthickness = 1; % fixed or variable layer thickness run

% Shelf forcing parameters
p.z0 = 50; % sets strength of shelf stratification
p.zd = 0; % sets strength of shelf oscillation, zero if no oscillation in shelf
p.tw = 0; % oscillation period of shelf forcing (days)
p.sf = @(S1, S2, Z, Z0) S1 - (S1 - S2)*exp(Z./Z0); % functional form of the shelf stratification
p.Stop = 30; % shelf salinity at the ocean surface
p.Sbottom = 35; % shelf salinity at the ocean floor
p.Ttop = 2.5; % shelf temperature at the ocean surface
p.Tbottom = 2.5; % shelf temperature at the ocean floor

% Glacier forcing parameters
p.Qv0 = 0; % volumetric flow rate of subglacial discharge
p.P0 = 25; % plume entrainment width, 0 = no plume
p.gf = @(T, Q) Q*ones(size(T)); % functional form of subglacial discharge rate depending on time

% Iceberg forcing parameters
p.icestatic = 1; % 1 for static icebergs, 0 for evolving
p.M0 = 7e-7; % iceberg melt efficiency
p.A0 = 2e9; % scaling of iceberg area
p.U0 = 0.25; % scale upwelling
p.nu0 = 25; % iceberg profile coefficient
% p.E0 = 1e-7; % iceberg export efficiency
% p.uIce = 0.005; % iceberg down-fjord velocity
p.if = @(NU, H, Z) (NU/H)*exp(NU*Z/H)/(1-exp(-NU)); % functional form of iceberg depth profile

%% Set the run output structure containing the inputted parameters.
run_output.p = p;
run_output.t = t;
% Set the boundary conditions for the input parameters at each timestep.
run_output.f = get_idealised_forcing(p, t);
% Set the initial conditions for the input parameters with the initial
% boundary conditions.
run_output.a = get_initial_conditions(p, run_output.f);

end