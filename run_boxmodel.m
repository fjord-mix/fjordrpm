%% run_boxmodel.m
% Description here

%%
clearvars
close all

%% Name the experiment and create an output directory

%% Load the problem setup file (parameters, initial conditions, forcings)
%% Examples
addpath('examples/');
% Example 1: intermediate_circulation
% Example 2:
% Example 3:
p = intermediate_circulation;

%% time steps etc
t = 0:p.dt:p.t_end;
f = get_forcing(p, t);
a = get_initial_conditions(p, f);

%% Run the model
%     % run model
% p.plot_runtime = 0;
% s = boxmodel_v4(p,f,a,t);
% 
% 
% % animate
% 
% %% and save the output
% save INT.mat s p f
% 
% %% visualisations
% 
% animate_v4p1('INT.mat','INT_N6_withsill',50);