% RUN_ZMODEL_PARALLEL Driver to run the zmodel simulation in parallel with
% default parameters example varying the magnitude of the subglacial
% discharge.

%% Dealing with paths
addpath(genpath('./'))
clearvars
close all
% Choose where to save your model outputs.
output_folder='./outputs/model_results/';

%% Set how many processes in parallel will be run.
num_workers=2;

%% Setup the model run.
% Get the model parameters, initial conditions and boundary conditions for
% the chosen example.
[p0, t, f, a] = get_model_default_parameters;
name = 'example_zmodel_parallel';
% Choose the parameter to vary and the values it will take (code will run
% in parallel for different values of this parameter).
parameter_to_vary = 'Qv0';
parameter_space = linspace(350, 1000, 50);

% Create an empty output structure to store the results.
if exist('fjord_par_outputs','Var'), clear fjord_par_outputs; end
fjord_par_outputs(length(parameter_space)) = struct("p",[],"a",[],"f",[],"t",[],"m",[],"s",[]);

% Set up the runs in serial to save memory.
for INDEX = 1:length(parameter_space)
    fjord_par_outputs(INDEX).m.name = sprintf('Iteration_%d\n', INDEX);
    fjord_par_outputs(INDEX).t = t;
    fjord_par_outputs(INDEX).p = mod_run_param(parameter_to_vary, parameter_space(INDEX), p0);
    % Compute new boundary conditions and initial conditions based on the
    % updated parameters p. Boundary conditions:
    fjord_par_outputs(INDEX).f = get_idealised_forcing(fjord_par_outputs(INDEX).p, t);
    % Initial conditions:
    fjord_par_outputs(INDEX).a = get_initial_conditions(fjord_par_outputs(INDEX).p, fjord_par_outputs(INDEX).f);
end

%% Run the model (output is saved automatically).
% If there is already a parallel pool running, do not create one.
checkPool = gcp('nocreate');
% If there is no pool, create one with the specified number of workers.
if isempty(checkPool)
    parpool(num_workers);
end

% Run the zmodel in parallel.
parfor INDEX = 1:length(parameter_space)
    fjord_par_outputs(INDEX).s = ...
        zmodel(fjord_par_outputs(INDEX).p, fjord_par_outputs(INDEX).t,fjord_par_outputs(INDEX).f,fjord_par_outputs(INDEX).a,[output_folder, fjord_par_outputs(INDEX).m.name, '.mat']);
    disp(fjord_par_outputs(INDEX).m.name)
end
