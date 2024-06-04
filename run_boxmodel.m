% RUN_BOXMODEL Driver to run the box model simulation.

%% Dealing with paths
addpath(genpath('./'))
clearvars
close all
output_folder='./outputs'; % choose where to save your model outputs here


%% Choosing the model example
% Setup for the model run. Load the user-defined run parameters. 
% Examples:
% 1. Intermediate circulation
% 2. Nudging based on target layers based on layer salinity
% 3. Subglacial discharge
% 4. Data-driven example: Kangerlussuaq during 2010-2018
% 5. Data-driven example: Benchmark fjords during 2010-2018 using 3 layers
%    above the sill (if exists).
% 6. Data-driven example: Benchmark fjords during 2010-2018 using 2 layers
%    above the sill (if exists).
% 7. Data-driven example: same as 5, but not taking care of some model
%    parametres (e.g., number of layers), which eventually causes the model
%    to crash. Not everything is perfect!
%
%    The last three examples have 4 fjords ready to be run: 
%    (1) Kangerlussuaq, (2) Sermilik, (3) Kangersuneq, (4) Ilulissat
%    Be mindful that these are not 100% "plug and play". The user
%    needs to choose which fjord from the fjord_model array to run

example_run = 5;
which_fjord = 3; % used for example_run 5 to 7

switch example_run
    case 1
        fjord_run = example_intermediate_circulation;
        name = 'example_intermediate_circulation';
    case 2
        fjord_run = example_nudging;
        name = 'example_nudging';
    case 3
        fjord_run = example_subglacial_discharge;
        name = 'example_subglacial_discharge';
    case 4
        fjord_run = load('./input_data_examples/KF_ctrl.mat').fjord_ctrl;
        name = 'example_Kangerlussuaq_2010_2018_3layers';
    case 5
        fjord_run = load('./input_data_examples/example_benchmark_fjords_3layers.mat').fjord_model(which_fjord);
        name = ['example_',fjord_run.m.name,'_2010_2018_3layers'];
    case 6        
        fjord_run = load('./input_data_examples/example_benchmark_fjords_2layers.mat').fjord_model(which_fjord);
        name = ['example_',fjord_run.m.name,'_2010_2018_2layers'];
    case 7
        fjord_run = load('./input_data_examples/example_benchmark_fjords_bad_3layers.mat').fjords_bad(which_fjord);
        name = ['bad_',fjord_run.m.name];
    otherwise
        fjord_run = test_changes;
        name = 'test_fjord_model';
end

% use the plot runtime at your discretion. It substantially slows down the simulation, because it spends most of the time plotting!
fjord_run.p.plot_runtime = 1; 

fjord_run.m.name = name;
[fjord_run.s,fjord_run.f] = boxmodel(fjord_run.p, fjord_run.t, fjord_run.f, fjord_run.a);

%% Saving results
if not(isfolder(output_folder))
    mkdir(output_folder)
    mkdir([output_folder,'/model_results']);
    mkdir([output_folder,'/figures']);
    mkdir([output_folder,'/animations']);
end

% Save the fjord structure, including input parameters, initial conditions, and results
save([output_folder,'/model_results/', name, '.mat'], 'fjord_run')

% create animation of the fjord run
animate([output_folder,'/model_results/', name, '.mat'],[output_folder,'/animations/'],name,50)

% Example plot that can be generated from the model outputs
plot_outputs(fjord_run);
exportgraphics(gcf,[output_folder,'/figures/',name,'.pdf'],'ContentType','vector','BackgroundColor','none')

plot_TS_boxmodel(fjord_run);
exportgraphics(gcf,[output_folder,'/figures/TS_',name,'.pdf'],'ContentType','vector','BackgroundColor','none')