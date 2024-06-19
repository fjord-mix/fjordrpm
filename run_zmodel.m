% RUN_ZMODEL Driver to run the zmodel simulation.

%% Dealing with paths
addpath(genpath('./'))
clearvars
close all
% Choose where to save your model outputs
output_folder='./outputs';

%% Choosing the model example
% Setup for the model run. Load the user-defined run parameters.
% Examples:
% 1. Intermediate circulation
% 2. Icebergs
% 3. Subglacial discharge
% 4. Data-driven example: something
%
% The data-driven example has 4 fjords ready to be run:
% (1) Kangerlussuaq, (2) Sermilik, (3) Kangersuneq, (4) Ilulissat
% The user needs to choose which fjord from the fjord_model array to run.

example_run = 10;
which_fjord = 3; % used for example_run 4

switch example_run
    case 1
        fjord_run = example_intermediate_circulation;
        name = 'example_intermediate_circulation';
    case 2
        fjord_run = example_icebergs;
        name = 'example_icebergs';
    case 3
        fjord_run = example_subglacial_discharge;
        name = 'example_subglacial_discharge';
    case 4
        fjord_run = load('a_file').fjord_ctrl;
        name = 'data_example';
    case 10 % debugging
        load('boxmodel_example_bad_H_negT_lim.mat');
        fjord_run = cur_fjord;
        name = 'martim_crash';
    otherwise
        fjord_run = test_changes;
        name = 'test_fjord_model';
end

%% Running the zmodel
% Use the plot runtime at your discretion. It substantially slows down the
% simulation, because it spends most of the time plotting!
fjord_run.p.plot_runtime = 0;
fjord_run.m.name = name;
[fjord_run.s,fjord_run.f] = zmodel(fjord_run.p, fjord_run.t, fjord_run.f, fjord_run.a);

%% Saving results
if not(isfolder(output_folder))
    mkdir(output_folder)
    mkdir([output_folder,'/model_results']);
    mkdir([output_folder,'/figures']);
    mkdir([output_folder,'/animations']);
end

% Save the fjord structure, including input parameters, initial conditions,
% and results
save([output_folder,'/model_results/', name, '.mat'], 'fjord_run')

%% Plotting results
% Create animation of the fjord run
% animate([output_folder,'/model_results/', name, '.mat'],[output_folder,'/animations/'],name,50)

% Example plots that can be generated from the model outputs
% plot_outputs(fjord_run);
% exportgraphics(gcf,[output_folder,'/figures/',name,'.pdf'],'ContentType','vector','BackgroundColor','none')
