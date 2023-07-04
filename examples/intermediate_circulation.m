%% intermediate_circulation.m
% Description here

%% Physical constants (these should not be changed in general)
p.g = 9.81; % gravity
p.betaS = 7.86e-4; % haline contraction
p.betaT = 3.87e-5; % thermal expansion
p.l = 3.34e5; % latent heat
p.cw = 3974; % water heat capacity
p.l1 = -5.73e-2; % usual linear freezing point constants
p.l2 = 9.32e-2;
p.l3 = -7.53e-4;
p.sid = 86400; % useful to have seconds in a day

%% Fjord geometry
p.W = 6e3; % fjord width
p.L = 80e3; % fjord length
p.H = 800; % fjord depth
p.zgl = -800; % grounding line depth
p.N = 6; % number of above-sill model layers
p.sill = 0; % flag for sill (0 = no sill, 1 = sill)
p.silldepth = 500; % only used if p.sill=1

%% Controlling and artificial parameters
p.P0 = 0; % no plume
p.C0 = 1e4;
p.K0 = 0.05;
p.M0 = 0; % no icebergs
p.nu0 = 25;
p.E0 = 1e-7;
p.dt = 0.1; % time stepping (units are days)
p.tmax = 1000;
t = 0:p.dt:p.tmax;

p.trelax = NaN; % no nudging
p.wmax = 4e-5;
p.Hmin = NaN; % no min thickness

%% Initial fjord conditions- this should be a separate function- get a? 
    % initial thicknesses
    if p.sill
        a.H0 = [(p.silldepth/p.N)*ones(1,p.N),p.H-p.silldepth];
    else
        a.H0 = [(p.H/p.N)*ones(1,p.N)];
    end

        % initial icebergs
    f.zi = [-800:10:0]';
    a.I0 = 0*f.zi;

    % forcings in structure f
    % icebergs
    f.D = zeros(1,length(t));
    f.xi = (p.nu0/sum(a.H0))*exp(p.nu0*f.zi/sum(a.H0))/(1-exp(-p.nu0));
    % subglacial discharge
    f.Qsg = zeros(1,length(t));
%% Initial shelf conditions- this should be a separate function- get f? 
% shelf stratification - oscillating
    f.zs = [-p.H:0];
    Sbottom = 35;
    Stop = 30;
    tw = 10; % 10 day period
    z0 = 50;
    zd = 30;
    zi = z0+(zd/2)*sin(2*pi*t/tw);
    for j=1:length(t),
        f.Ss(:,j) = Sbottom-(Sbottom-Stop)*exp(f.zs/zi(j));
        f.Ts(:,j) = 0*f.zs+3;
    end

        % set initial fjord T/S to be in equilibrium with shelf
    ints = [0;-cumsum(a.H0')];
    zs0 = unique(sort([f.zs,-cumsum(a.H0)]));
    Ss0 = interp1(f.zs,f.Ss(:,1),zs0,'pchip','extrap');
    Ts0 = interp1(f.zs,f.Ts(:,1),zs0,'pchip','extrap');
    for k=1:length(ints)-1,
        inds = find(zs0<=ints(k) & zs0>=ints(k+1));
        a.S0(k) = trapz(zs0(inds),Ss0(inds))/a.H0(k);
        a.T0(k) = trapz(zs0(inds),Ts0(inds))/a.H0(k);
    end