% PLOT_BOXMODEL Driver to plot the outputs of a particular box model simulation.
%% Setup for the plots
addpath('plotting_code/');
% Select the box model simulation.
name = 'EX_SGD';

%% 1. Box model animation
% Create a video of the box model simulation showing temperature and salinity
% and save to the output directory for the experiment.
nframes = 50; % number of frames in the animation
animate(name, nframes)


