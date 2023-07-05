% RUN_BOXMODEL Driver to run the box model simulation.
clearvars
close all

%% Get the run parameters
% Examples:
addpath('examples/');
% Example 1: intermediate_circulation
% Example 2:
% Example 3:
p = example_intermediate_circulation;

% Name the experiment and create a directory to store input parameters and
% outputs
name = 'EX_INT';
mkdir(['./output_',name]);

% Get the time variable, forcing structure and initial conditions structure
t = 0:p.dt:p.t_end;
f = get_forcing(p, t);
a = get_initial_conditions(p, f);

% Save the input parameters
save(['./output_',name,'/run_params.mat'], 'p', 'f', 't')

%% Run the model (output is saved automatically)
p.plot_runtime = 0;
% loop through INDEX to allow multiple box model runs (i.e. to study parameter space)
% can be a parfor loop
for INDEX = 1:1 
    boxmodel(p, f, a, t, name, INDEX);
end
