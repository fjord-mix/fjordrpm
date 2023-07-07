% RUN_BOXMODEL Driver to run the box model simulation.
clearvars
close all
addpath('simulation_code/')

% Setup for the model run. Load the user-defined run parameters. 
% Examples:
addpath('example_run_parameters/');
% 1. example_intermediate_circulation
% 2. example_nudging
% 3. example_subglacial_discharge 
%[p, a, f, t] = example_subglacial_discharge;
% 4. example_data_driven
% addpath('input_data_folder');
% my_nice_data = process_data(input_data);
% [p, a, ,f t] = example_data_driven(my_nice_data);

% Tests:
addpath('test_run_parameters/')
[p, a, f, t] = test_changes;
% Name the experiment and create a directory to store input parameters and
% outputs.
name = 'TEST_CHANGES';
mkdir(['./output_', name]);

% Save the input parameters and initial conditions
save(['./output_', name, '/run_params.mat'], 'p', 'a')

% Run the model (output is saved automatically).
p.plot_runtime = 0;
boxmodel(p, f, a, t, name);
