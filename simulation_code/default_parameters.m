function p = default_parameters

% DEFAULT_PARAMETERS Loads default parameters for fjordRPM.
%   p = DEFAULT_PARAMETERS sets the default box model parameters. The 
%   values can be modified afterwards for specific runs.

%% Physical constants - these should not be changed in general

p.g = 9.81;        % gravity (ms^-2)
p.betaS = 7.86e-4; % haline contraction (ppt^-1)
p.betaT = 3.87e-5; % thermal expansion (degC^-1)
p.l = 3.34e5;      % latent heat (J kg^-1)
p.cw = 3974;       % water heat capacity (J kg^-1 degC^-1)
p.sid = 86400;     % seconds in a day (useful to have)
p.alphai = 0.1;    % iceberg plume entrainment coefficient
p.alphap = 0.1;    % discharge plume entrainment coefficient

% Usual linear freezing point constants
p.l1 = -5.73e-2;   % dependence of freezing point on salinity (degC ppt^-1)
p.l2 = 9.32e-2;    % freezing point offset (degC)
p.l3 = -7.53e-4;   % dependence of freezing point on depth (degC m^-1)

%% Controlling parameters

p.wp = 250;        % plume width (m)
p.C0 = 1e4;        % shelf exchange efficiency (s)
p.K0 = 5e-3;       % vertical mixing efficiency (-)
p.Ri0 = 700;       % Richardson number dependency of mixing
p.M0 = 5e-7;       % iceberg melt efficiency (s^-1 degC^-1)
p.U0 = 1;          % scale iceberg upwelling

%% Run-time plotting
p.plot_runtime = 0;

end