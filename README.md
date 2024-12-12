# FjordRPM model repository

This is the official repository for the Fjord Reduced Physics Model (FjordRPM). 
A very brief explanation of the model is provided below, but all details can be read in the model description paper: [Slater et al., 2025 (dummy link)](https://github.com/fjord-mix/fjordrpm/)

## Installation

Just download or clone the GitHub repository by using `git clone https://github.com/fjord-mix/fjordrpm.git`

## Directory structure

- [`simulation_code`](https://github.com/fjord-mix/fjordrpm/tree/readme_update/simulation_code) contains the source code for FjordRPM
- [`plotting_code`](https://github.com/fjord-mix/fjordrpm/tree/readme_update/plotting_code) contains basic routines for plotting the raw results.
- [`examples`](https://github.com/fjord-mix/fjordrpm/tree/readme_update/examples) is a folder containing example simulations for the model, and are good templates to start your own simulations from

## Model usage

The model can be run by calling the function `s = run_model(p,t,f,a)`, which takes the following arguments:

- `t` the time axis for the simulation
    - `t(1)` should have the first time (usually, but not necessarily zero)	
    - `t(end)` should be the last time
    - the model time step is calculated as `t(j+1) - t(j)`

- `p`: the structure containing all model parameters. 
	- Default values for all physical parameters are provided in `default_parameters`, e.g., `p = default_parameters();`
	- Model geometry variables also need to go in p, e.g.:
	```
	p.Hgl = 500;   % grounding line depth (m, positive)
	p.Hsill = 650; % sill depth (m, positive)
	p.H = 1000;    % fjord depth (m, positive)
	p.L = 90e3;    % fjord length (m)
	p.W = 3e3;     % fjord width
	```
	- **Note 1:** `p.sill==1` tells the model the fjord has a sill, and `p.sill==0` tells the model the fjord doesn't. However, if `p.Hsill >= p.H` the model will automatically adjust to no sill. If `p.sill==0`, `p.Hsill` will be ignored.
	- **Note 2:** the user can define `p.t_save` as a subset of the original time axis `t` to reduce the number of time steps saved. If not specified, the model will save all time steps.
	
- `f`: the structure containing all model forcings
    - `f.Ts` is the shelf temperature (in $^oC$)
    - `f.Ss` is the shelf salinity (in PSU)
    - `f.zs` is the depths (in m) at which `Ts` and `Ss` are defined
    - `f.ts` is the time axis (usually in days, just like `t`) at which `Ts` and `Ss` are defined
    - `f.Qsg` is the subglacial discharge (in $m^3 s^{-1}$)
    - `f.tsg` is the time axis at which `Qsg` is defined, following the same reasoning as `f.ts`
    - **Note 1:** `f.Ts` and `f.Ss` must have dimensions `[length(f.zs), length(f.ts)]`
    - **Note 2:** `f.Qsg` must have dimensions `[n_plumes, length(f.tsg)]`. If `n_plumes > 1`, then `p.Hgl` needs to be an array containing the different grounding-line depths for each plume, such as `n_plumes == length(p.Hgl)`
    - **Note 3:** the model will automatically interpolate the forcing variables into the model's simulation time axis through the provided inputs `f.ts` and `f.tsg`
    
- `a`: the structure containing the initial conditions
	- `a.T0` contains the initial temperature profile
	- `a.S0` contains the initial salinity profile
	- `a.I0` contains the iceberg area profile for each layer
	- **Note:** all structures should have dimensions `[length(f.zs),1]`
	
- `s` is the output structure containing all variables
	- the variable `s.z` contains the depth axis of the fjord itself
	- the variable `s.t` contains the time axis at the time steps specified in `p.t_save`
	- dimensions of all variables are `[length(s.z),length(s.t)]`
	
**Note:** if a directory string is provided such as `s = run_model(p,t,f,a,'directory_path/filename.mat')`, then all inputs and outputs will be saved in `directory_path/filename.mat`