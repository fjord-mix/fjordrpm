% Script to demonstrate a FjordRPM simulation of fjord response to
% seasonally-varying subglacial discharge input. The shelf conditions are 
% constant in time and depth and there are no icebergs. The main example
% here assumes a single glacier with a single plume but the comments also
% show how to have multiple glaciers or plumes.

% clear workspace and close any figures to ensure clean environment
clear; close all;

% put FjordRPM code on path - need to update to the location of your code
path2sourcecode = '~/OneDrive - University of Edinburgh/fjordMIX/code/box-model/';
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
% (only used if there is non-zero subglacial discharge)
p.Hgl = 800; % grounding line depth (m)
p.Wp = 250; % subglacial discharge plume width (m)

% in the case of multiple plumes, either at the same glacier or at
% different glaciers, specify vectors of grounding line depth and plume
% width. for example, 3 glaciers of grounding line depth 800, 700, 600 m
% and plume width 300, 200 and 300 m would require
% p.Hgl = [800,700,600];
% p.Wp = [300,200,300];

% set up model layers
p.N = 40; % number of layers
a.H0 = (p.H/p.N)*ones(p.N,1); % layer thicknesses, here taken to be equal

% set up time stepping
dt = 0.2; % time step (in days)
t_end = 3*365; % time to end the simulation (in days)
t = 0:dt:t_end; % resulting time vector for simulation
p.t_save = 0:1:t_end; % times on which to save output

% set up shelf forcing - here constant in time and depth
% for more complexity see other examples
% f.ts must have dimensions 1 x nt
% f.zs must have dimensions nz x 1
% f.Ss and f.Ts must have dimensions nz x nt
f.ts = [0,t_end]; % time vector for shelf forcing
f.zs = [-p.H;0]; % depth vector for shelf forcing (negative below surface)
f.Ss = 34*ones(length(f.zs),length(f.ts)); % shelf salinity on (zs,ts)
f.Ts = 3*ones(length(f.zs),length(f.ts)); % shelf temperature on (zs,ts)

% set up subglacial discharge forcing
% here use idealised seasonal gaussian peaked at julian day 200
% f.tsg must have dimensions 1 x nt
% f.Qsg must have dimensions num plumes x nt
f.tsg = t; % time vector for subglacial discharge
f.Qsg = 300*exp(-((mod(t,365)-200)/30).^2); % subglacial discharge on tsg

% in the case of multiple plumes we add more rows to f.Qsg
% for example for three plumes with the same subglacial discharge
% f.Qsg(1,:) = 300*exp(-((mod(t,365)-200)/30).^2);
% f.Qsg(2,:) = 300*exp(-((mod(t,365)-200)/30).^2);
% f.Qsg(3,:) = 300*exp(-((mod(t,365)-200)/30).^2);

% fjord initial conditions
% set up to be same as initial shelf profiles
[a.T0, a.S0] = bin_shelf_profiles(f.Ts(:,1), f.Ss(:,1), f.zs, a.H0);

% set up icebergs - in this example there are no icebergs
a.I0 = 0*a.H0;

% run the model
% p.plot_runtime = 1; % plot while simulation runs - fun but quite slow
s = run_model(p, t, f, a);

% save the output
save example_subglacial_discharge.mat s p t f a

% make an animation of the output (takes a few minutes)
% animate(p,s,50,'example_subglacial_discharge');

% make basic plots of the output
plotrpm(p,s);




