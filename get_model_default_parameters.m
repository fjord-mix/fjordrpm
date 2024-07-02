function [p, t, f, a] = get_model_default_parameters

% GET_MODEL_DEFAULT_PARAMETERS  Loads default parameters for zmodel.
%   [P, A] = GET_MODEL_DEFAULT_PARAMETERS sets the default box model
%   parameters for most examples. The values can be modified afterwards for
%   specific runs.

%% Physical constants - these should not be changed in general
p.g = 9.81;        % gravity (ms^-2)
p.betaS = 7.86e-4; % haline contraction (ppt^-1)
p.betaT = 3.87e-5; % thermal expansion (deg.C^-1)
p.l = 3.34e5;      % latent heat (J kg^-1)
p.cw = 3974;       % water heat capacity (J kg^-1 deg.C^-1)
p.sid = 86400;     % seconds in a day (useful to have)

% Usual linear freezing point constants
p.l1 = -5.73e-2;   % dependence of freezing point on salinity (deg.C ppt^-1)
p.l2 = 9.32e-2;    % freezing point offset (deg.C)
p.l3 = -7.53e-4;   % dependence of freezing point on depth (deg.C m^-1)

%% Controlling parameters
p.P0 = 25;         % entrainment efficiency (m)
p.C0 = 1e4;        % shelf exchange efficiency (s)
p.K0 = 5e-3;       % vertical mixing efficiency (-)
p.Ri0 = 700;       % Richardson number dependency of mixing
p.M0 = 2e-8;       % iceberg melt efficiency (s^-1 deg.C^-1)
p.nu0 = 25;        % iceberg volume profile coefficient (-)
p.E0 = 1e-7;       % iceberg export efficiency (s^-1)
p.uIce = 0.005;    % iceberg down-fjord velocity (m s^-1)
p.gamma = 0.5;     % proportion of iceberg melt flux that gets mixed vertically (-)
p.alphaI = 0.1;    % iceberg plume entrainment coefficient
p.icestatic = 1;   % whether icebergs are "static" or "dynamic" in the model
p.A0 = 2e9;        % scaling of iceberg area
p.U0 = 0.25;       % scale upwelling

%% Idealised fjord geometry
p.W = 6e3;          % fjord width
p.L = 80e3;         % fjord length
p.H = 800;          % fjord depth
p.silldepth = -500; % only used if p.sill=1
p.zgl = -800;       % grounding line depth
p.Hmin = 5;        % min box thickness

%% Idealised-forcing parameters
%  These should only be used for idealised boundary conditions.

% Shelf forcing parameters
p.z0 = 50; % sets strength of shelf stratification
p.zd = 30; % sets strength of shelf oscillation, zero if no oscillation in shelf
p.tw = 10; % oscillation period of shelf forcing (days)
p.Stop = 30; % shelf salinity at the ocean surface
p.Sbottom = 35; % shelf salinity at the ocean floor
p.Ttop = 3; % shelf temperature at the ocean surface
p.Tbottom = 3; % shelf temperature at the ocean floor
p.sf = @(S1, S2, Z, Z0) S1 - (S1 - S2)*exp(Z./Z0); % functional form of the shelf stratification

% Glacier forcing parameters
p.Qv0 = 0; % volumetric flow rate of subglacial discharge
p.gf = @(T, Q) Q*ones(size(T)); % functional form of subglacial discharge rate depending on time

% Iceberg forcing parameters
p.D0 = 0; % volumetric iceberg discharge from glacier(s) into the fjord
p.if = @(NU, H, Z) (NU/H)*exp(NU*Z/H)/(1-exp(-NU)); % functional form of iceberg depth profile

%% ZMODEL properties
p.N    = 60;    % 3 layers above the sill
p.sill = 1;    % and sill exists, totalling 4 layers

%% Parameters for plotting
p.plot_runtime=0;

%% Time values at which to compute the solution (in days)
t = 0:1:200;
% Time values at which to save the solution (in days)

%% Set the boundary and initial conditions
% Boundary conditions
f = get_idealised_forcing(p, t);
% Initial conditions
a = get_initial_conditions(p, f);

end