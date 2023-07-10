% RUN_BOXMODEL Driver to run the box model simulation.
addpath(genpath('./'))
clearvars
close all

% Setup for the model run. Load the user-defined run parameters. 
% Examples:
% 1. example_intermediate_circulation
% example_run = example_intermediate_circulation;
% 2. example_nudging
% example_run = example_nudging;
% 3. example_subglacial_discharge 
% example_run = example_subglacial_discharge;
% 4. example_data_driven
% addpath('input_data_folder');
% example_run = load('./example_input_data/KF_ctrl.mat').fjord_ctrl;

% Tests:
example_run = test_changes;
% Name the experiment and create a directory to store input parameters and
% outputs.
name = 'TEST_CHANGES';
mkdir(['./output_', name]);

% Save the input parameters and initial conditions
save(['./output_', name, '/run_params.mat'], 'p', 'a')

% Run the model (output is no longer saved automatically).
example_run.m.name=name;
p.plot_runtime = 1;
[example_run.s,example_run.f] = boxmodel(example_run.p, example_run.f, example_run.a, example_run.t);
