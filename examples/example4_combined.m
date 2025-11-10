% Script to demonstrate a 2-year FjordRPM simulation of a fjord with
% subglacial discharge, icebergs and shelf variability (i.e., something of
% a combination of the three examples that look at each individually). The
% subglacial discharge is seasonal, peaking at day 200 in each year. The
% icebergs are constant in time with a total iceberg-ocean surface area of
% 100 km^2. There is shelf variability only outside of summer.

% clear workspace and close any figures to ensure clean environment
clear; close all;

% put FjordRPM code on path - need to update to the location of your code
path2sourcecode = '~/OneDrive - University of Edinburgh/fjordMIX/code/fjordrpm/';
addpath(genpath(path2sourcecode));

% get basic constants and default controlling parameters
p = default_parameters;

% can adjust any of the default parameters afterwards if needed
% p.C0 = 5e4; % for example adjust shelf exchange parameter
% p.run_plume_every = 10; % or update plume model only every 10 time steps

% set up fjord geometry
p.W = 6e3; % fjord width (m)
p.L = 60e3; % fjord length (m)
p.H = 800; % fjord depth (m)
p.sill = 1; % p.sill=1 for presence of sill, p.sill=0 for no sill
p.Hsill = 400; % sill depth below surface (m), only used if p.sill=1

% set up glacier geometry
% (not actually used since there is no subglacial discharge)
p.Hgl = 800; % grounding line depth (m)
p.Wp = 250; % subglacial discharge plume width (m)

% set up model layers
p.N = 40; % number of layers
a.H0 = (p.H/p.N)*ones(p.N,1); % layer thicknesses, here taken to be equal

% set up time stepping
dt = 0.2; % time step (in days)
t_end = 2*365; % time to end the simulation (in days)
t = 0:dt:t_end; % resulting time vector for simulation
p.t_save = 0:1:t_end; % times on which to save output

% set up shelf forcing
% here make it an exponential profile that shoals and deepens every 10 days
% f.ts must have dimensions 1 x nt
% f.zs must have dimensions nz x 1
% f.Ss and f.Ts must have dimensions nz x nt
f.ts = 0:1:t_end; % time vector for shelf forcing
f.zs = [-p.H:0]'; % depth vector for shelf forcing (negative below surface)
Sbottom = 35; % salinity at bottom
Stop = 30; % salinity at top
Tbottom = 4; % temperature at bottom
Ttop = 0; % temperature at top
tw = 10; % period of oscillation (days)
zi = 50+(30/2)*sin(2*pi*f.ts/tw); % 'pycnocline' oscillation
zi(mod(f.ts,365)>120 & mod(f.ts,365)<280) = 50; % no oscillation in summer
for k=1:length(f.ts),
    f.Ss(:,k) = Sbottom-(Sbottom-Stop)*exp(f.zs/zi(k)); % shelf salinity
    f.Ts(:,k) = Tbottom-(Tbottom-Ttop)*exp(f.zs/zi(k)); % shelf temperature
end

% set up subglacial discharge forcing
% here use idealised seasonal gaussian peaked at julian day 200
% f.tsg must have dimensions 1 x nt
% f.Qsg must have dimensions num plumes x nt
f.tsg = t; % time vector for subglacial discharge
f.Qsg = 300*exp(-((mod(f.tsg,365)-200)/30).^2); % subglacial discharge on tsg

% fjord initial conditions
% set up to be same as initial shelf profiles
[a.T0, a.S0] = bin_shelf_profiles(f.Ts(:,1), f.Ss(:,1), f.zs, a.H0);

% set up icebergs
% need to specify the iceberg-ocean surface area in each layer. assume an 
% exponential profile with total surface area 100 km^2 = 1e8 m^2
zc = cumsum(a.H0)-a.H0/2; % depth of centre of layers
iceprofile = exp(-zc/100); % profile before scaling
a.I0 = (1e8/sum(iceprofile))*iceprofile; % surface area profile used by FjordRPM (m^2)

% run the model
% p.plot_runtime = 1; % plot while simulation runs - fun but quite slow
s = run_model(p, t, f, a);

% save the output
save example4_combined.mat s p t f a

% make an animation of the output (takes a few minutes)
% animate(p,s,100,'example4_combined');

% make basic plots of the output
plotrpm(p,s,50);





