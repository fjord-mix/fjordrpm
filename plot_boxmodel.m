% PLOT_BOXMODEL Driver to plot the outputs of a particular box model simulation.
clearvars
close all

%% Setup for the plots
addpath('plots/');
% Select the box model simulation.
name = 'EX_QSG';

%% 1. Box model animation
% Create a video of the box model simulation and save to the output
% directory for the experiment
nframeds = 50; % number of frames in the animation
animatep1(name, nframes)


