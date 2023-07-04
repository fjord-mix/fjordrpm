%% intermediate_circulation.m

function p = intermediate_circulation

% Description here
% p = parameters (user sets these)
% a = initial conditions
% f = forcing from shelf, icebergs and subglacial discharge

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
p.zd = 30; % sets strength of shelf oscillation, zero if no oscillation in shelf
p.tw = 10; % oscillation period of shelf forcing (days)
p.sf = @(S1, S2, Z, Z0) S1 - (S1 - S2)*exp(Z/Z0); % functional form of the shelf stratification

%% Glacier parameters
p.Qv0 = 0; % volumetric flow rate of subglacial discharge
p.P0 = 0; % plume entrainment width, 0 = no plume 

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

%% Time step
p.dt = 0.1; % time stepping (units are days)
p.t_end = 1000;





%% Forcing parameters (shelf, iceberg and subglacial discharge)
% shelf




% function f = get_forcing(p)
% shelf
% f.zs = -p.H:0; % will this break if p is not an integer?
% f.zi = p.z0+(p.zd/2)*sin(2*pi*t/p.tw);
% t = 0:p.dt:p.t_end;
% for j=1:length(t),
%     f.Ss(:,j) = p.sf(p.Sbottom, p.Stop, f.zs, f.zi(j));
%     f.Ts(:,j) = 0*f.zs+3;
% end
%icebergs
%f.zi = [-800:10:0]'; % should be in terms of paramaters
% icebergs
% f.D = zeros(1,length(t));
% f.xi = (p.nu0/sum(a.H0))*exp(p.nu0*f.zi/sum(a.H0))/(1-exp(-p.nu0));
% % subglacial discharge
% % subglacial discharge
% f.Qsg = zeros(1,length(t));
%  f.Qsg(find(t>5)) = 100;
% icebergs
%end


%% Initial fjord conditions
%function a = get_initial_conditions(p, f)
% initial thicknesses
if p.sill
    a.H0 = [(p.silldepth/p.N)*ones(1,p.N),p.H-p.silldepth];
else
    a.H0 = [(p.H/p.N)*ones(1,p.N)];
end

% initial icebergs

a.I0 = 0*f.zi;


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

end