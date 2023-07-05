function p = example_subglacial_discharge

% EXAMPLE_SUBGLACIAL_DISCHARGE  Example to run the box model with
% subglacial discharge.
%   P = EXAMPLE_SUBGLACIAL_DISCHARGE sets the user-defined box model parameters
%   P for the subglacial discharge example.

%% Physical constants (these should not be changed in general)
p.g = 9.81; % gravity
p.betaS = 7.86e-4; % haline contraction
p.betaT = 3.87e-5; % thermal expansion
p.l = 3.34e5; % latent heat
p.cw = 3974; % water heat capacity
p.l1 = -5.73e-2; % usual linear freezing point constants
p.l2 = 9.32e-2;
p.l3 = -7.53e-4;
p.sid = 86400; % seconds in a day

%% Fjord geometry parameters
p.W = 6e3; % fjord width
p.L = 80e3; % fjord length
p.H = 800; % fjord depth
p.zgl = -800; % grounding line depth
p.N = 6; % number of above-sill model layers
p.sill = 0; % flag for sill (0 = no sill, 1 = sill)
p.silldepth = 500; % only used if p.sill=1

%% Shelf parameters
p.Stop = 30; % shelf salinity at the ocean surface
p.Sbottom = 35; % shelf salinity at the ocean floor
p.z0 = 50; % sets strength of shelf stratification
p.zd = 0; % sets strength of shelf oscillation, zero if no oscillation in shelf
p.tw = 0; % oscillation period of shelf forcing (days)
p.sf = @(S1, S2, Z, Z0) S1 - (S1 - S2)*exp(Z./Z0); % functional form of the shelf stratification

%% Glacier parameters
p.Qv0 = 100; % volumetric flow rate of subglacial discharge
p.P0 = 25; % plume entrainment width, 0 = no plume

%% Iceberg parameters
p.M0 = 0; % iceberg melt efficiency, 0 = no icebergs
p.nu0 = 25; % iceberg volume profile coefficient
p.E0 = 1e-7; % iceberg export efficiency

%% Tuning/artificial parameters
p.C0 = 1e4; % shelf exchange effiency
p.K0 = 0.05; % vertical mixing effiency
p.wmax = 4e-5; % sets an upper bound on vertical mixing
p.Hmin = NaN; % minimum box thickness, NaN = no min thickness
p.trelax = NaN; % controls layer nudging, NaN = no nudging
p.Snudge = 0.*ones(1, p.N-1); % controls layer thickness for nudging, 0 if no nudging

%% Time step
p.dt = 0.1; % time stepping (units are days)
p.t_end = 10; % time to end the simulation

end