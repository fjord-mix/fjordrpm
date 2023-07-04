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
% 
%% and save the output 
