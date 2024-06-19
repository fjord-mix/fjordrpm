% RUN_ZMODEL Driver to run the zmodel simulation in parallel.
% Examples: SGD-driven, constant shelf temp, stratified shelf salinity.

num_workers=2; % how many processes in parallel will be run
%% Set up parameter space
experiment = 'sgd';
parameter_to_vary = 'Qv0';
parameter_space = linspace(350, 1000, 50);
t = linspace(1,100,1);

if exist('fjord_par_outputs','Var'), clear fjord_par_outputs; end
fjord_par_outputs(length(parameter_space)) = struct("p",[],"a",[],"f",[],"t",[],"m",[],"s",[]);

p_def = get_model_default_parameters; % Get default parameters
% setting up runs in serial significantly saves memory
for INDEX = 1:length(parameter_space)
    fjord_par_outputs(INDEX).m.name = sprintf('Iteration_%d\n', INDEX);
    fjord_par_outputs(INDEX).t = t;
    fjord_par_outputs(INDEX).p = mod_run_param(parameter_to_vary, parameter_space(INDEX), p_def);
    % Set idealised boundary and initial conditions, if not given, based on input parameters.
    % Boundary conditions:
    if isempty(fjord_par_outputs(INDEX).f)
        fjord_par_outputs(INDEX).f = get_idealised_forcing(fjord_par_outputs(INDEX).p, t);
    end % we cannot use an OR statement here
    % Initial conditions
    if isempty(fjord_par_outputs(INDEX).a)
        fjord_par_outputs(INDEX).a = get_initial_conditions(fjord_par_outputs(INDEX).p, fjord_par_outputs(INDEX).f);
    end
end

%% Run the model (output is saved automatically).
% starts parallel pool
checkPool = gcp('nocreate'); % If no pool, do not create one
if isempty(checkPool) % if there is no pool
    parpool(num_workers);
end

parfor INDEX = 1:length(parameter_space)
    [fjord_par_outputs(INDEX).s,fjord_par_outputs(INDEX).f] = ...
        zmodel(fjord_par_outputs(INDEX).p, fjord_par_outputs(INDEX).t,[],[],[outs_path,'test_parallel/',fjord_par_outputs(INDEX).m.name,'.mat']);
    disp(fjord_par_outputs(INDEX).m.name)
end