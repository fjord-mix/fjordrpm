% RUN_BOXMODEL Driver to run the box model simulation.
addpath(genpath('./'))
clearvars
close all

% Setup for the model run. Load the user-defined run parameters. 
% Examples:
% 1. example_intermediate_circulation
% 2. example_nudging
% 3. example_subglacial_discharge 
% run_subglacial_discharge = example_subglacial_discharge;
% 4. example_data_driven
% addpath('input_data_folder');
% kangerlussuaq_example = load('./example_input_data/KF_ctrl.mat').fjord_ctrl;
% kangerlussuaq_example.s = boxmodel(kangerlussuaq_example.p,kangerlussuaq_example.f,kangerlussuaq_example.a,kangerlussuaq_example.t);

% Tests:
test_fjord_run = test_changes;
% Name the experiment and create a directory to store input parameters and
% outputs.
name = 'TEST_CHANGES';
mkdir(['./output_', name]);

% Save the input parameters and initial conditions
save(['./output_', name, '/run_params.mat'], 'p', 'a')

% Run the model (output is saved automatically).
p.plot_runtime = 0;
boxmodel(p, f, a, t, name);
