% RUN_ZMODEL Driver to run the zmodel simulation with default model
% parameters, and save and plot the results.
addpath(genpath('./'))
clearvars
close all

% Get default model inputs.
[fjord_run.p, fjord_run.t, fjord_run.f, fjord_run.a] = get_model_default_parameters;
name = 'default_example';

% Choose whether to plot the simulation as it is running.
fjord_run.p.plot_runtime = 0;

% Run the ZMODEL.
fjord_run.s = zmodel(fjord_run.p, fjord_run.t, fjord_run.f, fjord_run.a);

% Choose where to save the model outputs.
output_folder='./outputs';
% Save the results.
if not(isfolder(output_folder))
    mkdir(output_folder)
    mkdir([output_folder,'/model_results']);
    mkdir([output_folder,'/figures']);
    mkdir([output_folder,'/animations']);
end

% Save the fjord structure, including input parameters, initial conditions,
% boundary conditions and results
save([output_folder, '/model_results', name, '.mat'], 'fjord_run')

% Plot the results.
plot_outputs(fjord_run);
exportgraphics(gcf,[output_folder,'/figures/',name,'.pdf'],'ContentType','vector','BackgroundColor','none')

% Create an animation of the fjord run.
% animate([output_folder,'/model_results/', name, '.mat'],[output_folder,'/animations/'],name,50)


