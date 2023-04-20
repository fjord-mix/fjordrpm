clear; close all;

% 1 for INT
% 2 for MINTHICKNESS
% 3 for LOW SUBGLACIAL DISCHARGE
% 4 for HIGH SUBGLACIAL DISCHARGE
% 5 for NUDGING
% 6 for NUDGING WITH HIGH SUBGLACIAL DISCHARGE
runno = 1;

% physical constants
p.g = 9.81; % gravity
p.betaS = 7.86e-4; % haline contraction
p.betaT = 3.87e-5; % thermal expansion
p.l = 3.34e5; % latent heat
p.cw = 3974; % water heat capacity
p.l1 = -5.73e-2; % usual linear freezing point constants
p.l2 = 9.32e-2;
p.l3 = -7.53e-4;
p.sid = 86400; % useful to have seconds in a day

% controlling parameters
p.P0 = 25;
p.C0 = 1e4;
p.K0 = 0.05;
p.M0 = 2e-8;
p.nu0 = 25;
p.E0 = 1e-7;

% artificial parameters
p.wmax = 4e-5;
p.S12 = 33;
p.S23 = 33.5;
p.trelax = 10*86400; % in seconds

% fjord geometry
p.W = 6e3; % fjord width
p.L = 80e3; % fjord length
p.H = 800; % fjord depth
p.zgl = -800; % grounding line depth
p.Hmin = 25; % min box thickness

%% shelf-driven (intermediary) circulation
if runno==1,

    % turn off processes not needed
    p.P0 = 0; % no plume
    p.M0 = 0; % no icebergs
    p.trelax = NaN; % no nudging

    % time stepping (units are days)
    p.dt = 0.1;
    t = [0:p.dt:100];

    % initial thicknesses
    a.H0 = [50,100,350,300];

    % initial icebergs
    f.zi = [-800:10:0]';
    a.I0 = 0*f.zi;

    % forcings in structure f
    % icebergs
    f.D = zeros(1,length(t));
    f.xi = (p.nu0/sum(a.H0))*exp(p.nu0*f.zi/sum(a.H0))/(1-exp(-p.nu0));
    % subglacial discharge
    f.Qsg = zeros(1,length(t));
    % shelf stratification - oscillating
    f.zs = [-500:0];
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

    % run model
    s = boxmodel_v4(p,f,a,t);
    save INT.mat s p f

    % animate
    animate_v4('INT.mat','INT',50);

end

%% min thickness illustration
if runno==2,

    % turn off processes not needed
    p.P0 = 0; % no plume
    p.M0 = 0; % no icebergs
    p.trelax = NaN; % no nudging
    
    % turn up shelf coefficient to try to make a layer disappear
    p.C0 = 1e5;

    % time stepping (units are days)
    p.dt = 0.1;
    t = [0:p.dt:100];

    % initial thicknesses
    a.H0 = [50,100,350,300];

    % initial icebergs
    f.zi = [-800:10:0]';
    a.I0 = 0*f.zi;

    % forcings in structure f
    % icebergs
    f.D = zeros(1,length(t));
    f.xi = (p.nu0/sum(a.H0))*exp(p.nu0*f.zi/sum(a.H0))/(1-exp(-p.nu0));
    % subglacial discharge
    f.Qsg = zeros(1,length(t));
    % shelf stratification - oscillating
    f.zs = [-500:0];
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

    % run model
    s = boxmodel_v4(p,f,a,t);
    save MINTHICKNESS.mat s p f

    % animate
    animate_v4('MINTHICKNESS.mat','MINTHICKNESS',50);

end

%% LOWQSG
if runno==3,

    % turn off processes not needed
    p.M0 = 0; % no icebergs
    p.trelax = NaN;

    % time stepping (units are days)
    p.dt = 0.25;
    t = [0:p.dt:800];

    % initial thicknesses
    a.H0 = [50,100,350,300];

    % initial icebergs
    f.zi = [-800:10:0]';
    a.I0 = 0*f.zi;

    % forcings in structure f
    % icebergs
    f.D = zeros(1,length(t));
    f.xi = (p.nu0/sum(a.H0))*exp(p.nu0*f.zi/sum(a.H0))/(1-exp(-p.nu0));
    % subglacial discharge
    f.Qsg = zeros(1,length(t));
    f.Qsg(find(t>10)) = 100;
    % shelf stratification - constant
    f.zs = [-500:0];
    Sbottom = 35;
    Stop = 30;
    z0 = 50;
    for j=1:length(t),
        f.Ss(:,j) = Sbottom-(Sbottom-Stop)*exp(f.zs/z0);
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

    % run model
    s = boxmodel_v4(p,f,a,t);
    save LOWQSG.mat s p f

    % animate
    animate_v4('LOWQSG.mat','LOWQSG',50);

end

%% HIGHQSG
if runno==4,

    % turn off processes not needed
    p.M0 = 0; % no icebergs
    p.trelax = NaN;

    % time stepping (units are days)
    p.dt = 0.25;
    t = [0:p.dt:800];

    % initial thicknesses
    a.H0 = [50,100,350,300];

    % initial icebergs
    f.zi = [-800:10:0]';
    a.I0 = 0*f.zi;

    % forcings in structure f
    % icebergs
    f.D = zeros(1,length(t));
    f.xi = (p.nu0/sum(a.H0))*exp(p.nu0*f.zi/sum(a.H0))/(1-exp(-p.nu0));
    % subglacial discharge
    f.Qsg = zeros(1,length(t));
    f.Qsg(find(t>10)) = 500;
    % shelf stratification - constant
    f.zs = [-500:0];
    Sbottom = 35;
    Stop = 30;
    z0 = 50;
    for j=1:length(t),
        f.Ss(:,j) = Sbottom-(Sbottom-Stop)*exp(f.zs/z0);
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

    % run model
    s = boxmodel_v4(p,f,a,t);
    save HIGHQSG.mat s p f

    % animate
    animate_v4('HIGHQSG.mat','HIGHQSG',50);

end

%% NUDGING
if runno==5,

    % turn off processes not needed
    P.P0 = 0; % no plume
    p.M0 = 0; % no icebergs

    % nudging
    p.trelax = 10*86400; % in seconds
    p.S12 = 34.3;
    p.S23 = 34.9;

    % time stepping (units are days)
    p.dt = 0.25;
    t = [0:p.dt:800];

    % initial thicknesses
    a.H0 = [50,100,350,300];

    % initial icebergs
    f.zi = [-800:10:0]';
    a.I0 = 0*f.zi;

    % forcings in structure f
    % icebergs
    f.D = zeros(1,length(t));
    f.xi = (p.nu0/sum(a.H0))*exp(p.nu0*f.zi/sum(a.H0))/(1-exp(-p.nu0));
    % subglacial discharge
    f.Qsg = zeros(1,length(t));
    % shelf stratification - constant
    f.zs = [-500:0];
    Sbottom = 35;
    Stop = 30;
    z0 = 50;
    for j=1:length(t),
        f.Ss(:,j) = Sbottom-(Sbottom-Stop)*exp(f.zs/z0);
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

    % run model
    s = boxmodel_v4(p,f,a,t);
    save NUDGING.mat s p f

    % animate
    animate_v4('NUDGING.mat','NUDGING',50);

end

%% NUDGING_HIGHQSG
if runno==6,

    % turn off processes not needed
    p.M0 = 0; % no icebergs

    % nudging
    p.trelax = 10*86400; % in seconds
    p.S12 = 34.3;
    p.S23 = 34.9;

    % time stepping (units are days)
    p.dt = 0.25;
    t = [0:p.dt:800];

    % initial thicknesses
    a.H0 = [50,100,350,300];

    % initial icebergs
    f.zi = [-800:10:0]';
    a.I0 = 0*f.zi;

    % forcings in structure f
    % icebergs
    f.D = zeros(1,length(t));
    f.xi = (p.nu0/sum(a.H0))*exp(p.nu0*f.zi/sum(a.H0))/(1-exp(-p.nu0));
    % subglacial discharge
    f.Qsg = zeros(1,length(t));
    f.Qsg(find(t>10)) = 500;
    % shelf stratification - constant
    f.zs = [-500:0];
    Sbottom = 35;
    Stop = 30;
    z0 = 50;
    for j=1:length(t),
        f.Ss(:,j) = Sbottom-(Sbottom-Stop)*exp(f.zs/z0);
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

    % run model
    s = boxmodel_v4(p,f,a,t);
    save NUDGING_HIGHQSG.mat s p f

    % animate
    animate_v4('NUDGING_HIGHQSG.mat','NUDGING_HIGHQSG',50);

end