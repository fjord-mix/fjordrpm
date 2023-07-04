%% run_boxmodel.m
% Description here 

%%
clearvars
close all

%% Name the experiment and create an output directory

%% Load the problem setup file (parameters, initial conditions, forcings)
addpath('../examples');
p = intermediate_circulation;
t = 0:p.dt:p.t_end;
f = get_forcing(p, t);
a = get_initial_conditions(p, f);

%% Run the model and save the output 
