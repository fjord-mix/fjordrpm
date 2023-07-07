function [p,a] = get_model_default_parameters()

% GET_MODEL_DEFAULT_PARAMETERS  Loads default parameters for boxmodel
%   [p,a] = GET_MODEL_DEFAULT_PARAMETERS sets the default box model parameters
%   for most examples. The values can be modified afterwards for specific runs

    %% Physical constants - these should not me changed in general
    p.g = 9.81;        % gravity (ms^-2)
    p.betaS = 7.86e-4; % haline contraction (ppt^-1)
    p.betaT = 3.87e-5; % thermal expansion (deg.C^-1)
    p.l = 3.34e5;      % latent heat (J kg^-1)
    p.cw = 3974;       % water heat capacity (J kg^-1 deg.C^-1)
    p.sid = 86400;     % seconds in a day (useful to have)
    
    % usual linear freezing point constants
    p.l1 = -5.73e-2;   % dependence of freezing point on salinity (deg.C ppt^-1)
    p.l2 = 9.32e-2;    % freezing point offset (deg.C)
    p.l3 = -7.53e-4;   % dependence of freezing point on depth (deg.C m^-1)    

    % controlling parameters
    p.P0 = 25;         % entrainment efficiency (m)
    p.C0 = 1e4;        % shelf exchange efficiency (s)
    p.K0 = 0.05;       % vertical mixing efficiency (-)
    p.M0 = 2e-8;       % iceberg melt efficiency (s^-1 deg.C^-1)
    p.nu0 = 25;        % iceberg volume profile coefficient (-)
    p.E0 = 1e-7;       % iceberg export efficiency (s^-1)

    % artificial parameters
    p.wmax = 4e-5;         % maximum vertical mixing velocity    
    p.trelax = NaN;        % nudging relaxation time in days
    p.real_time_nudge = 0; % set to 1 to update nudging values from current shelf conditions

    % density-based boundaries to define surface, intermediate, and deep layers in sigma(theta), kg m^-3
    % p.sigma_bnds = [27.3, 27.7]; % upper bnd according to Cowton et al., lower bnd with T > 0 from Rudels et al. (2002)
    p.sigma_bnds = [26.6, 27.3, 27.5, 27.6]; % adds dummy values in case of extra layers
    

    % idealised fjord geometry
    p.W = 6e3;          % fjord width
    p.L = 80e3;         % fjord length
    p.H = 800;          % fjord depth
    p.silldepth = -500; % only used if p.sill=1
    p.zgl = -800;       % grounding line depth
    p.Hmin = 25;        % min box thickness

    % initial thicknesses    
    p.N    = 3;    % 3 layers above the sill
    p.sill = 1;    % and sill exists, totalling 4 layers
    a.H0   = [50,100,350,300]; % thicknesses of all 4 layers

    p.plot_runtime=0;
    p.debug=0;
end