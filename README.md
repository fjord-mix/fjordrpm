
# Fjord Reduced Physics Model (FjordRPM)

The Fjord Reduced Physics Model (FjordRPM) is an efficient "1.5-dimensional" model for simulating the dynamics of glacial fjords. The model splits a fjord into a number of vertically-stacked layers and solves for the evolution of layer properties due to processes that currently include (i) plumes driven by subglacial discharge, (ii) exchange with the continental shelf, (iii) iceberg melt and upwelling, (iv) vertical mixing and advection within the fjord. The model can represent multiple plumes and also the presence of a sill. Full details on the physics, numerical implementation and validation of the model can be found in the model description paper. 

Here we provide instructions for downloading and running the model, together with a brief overview of the inputs/outputs and the provided examples. The model is written in Matlab and was tested mostly in R2022b but we do not expect difficulties with different versions of Matlab.

## Installation

Just download directly from GitHub/Zenodo or clone the GitHub repository by using `git clone https://github.com/fjord-mix/fjordrpm.git`. Ensure the code is on your Matlab path by executing `addpath(genpath(path2sourcecode))` in Matlab, where `path2sourcecode` is the location of the code on your machine. This is also done at the start of the code for the examples below but you will need to update for the location of the code on your machine.

## Directory structure

- `simulation_code` contains the source code for FjordRPM.
- `plotting_code` contains scripts giving examples of plotting and animating the results.
- `examples` contains example simulations that are a good place to start and to adapt to your own simulations.

## Model usage

The model is run by executing (in Matlab) `s = run_model(p,t,f,a)`, in which the arguments are

- `p` is a structure containing all the model parameters. 
    - Default values for most of the parameters are provided in `simulation_code/default_parameters` and can be loaded by executing `p = default_parameters()`
    - `p.N` defines the number of model layers.
    - The geometry of the glacier-fjord system is also specified in this structure, for example:
    ```
    p.Hgl = 500;   % glacier grounding line depth (m, positive)
    p.Hsill = 650; % sill depth (m, positive)
    p.H = 1000;    % fjord depth (m, positive)
    p.L = 90e3;    % fjord length (m)
    p.W = 3e3;     % fjord width (m)
    ```
    - Setting `p.sill==1` is needed to tell the model that the fjord has a sill. If we define `p.sill==0` then the model will assume there is no sill (even if `p.Hsill` is defined).
    - `p.run_plume_every` defines how often the plume model sitting inside FjordRPM is run. For example, `p.run_plume_every=20` means that the subglacial discharge plume dynamics are updated every 20 model time steps. Since the plume model is the most time-consuming aspect of the code, this can speed up the model.
    - The user can define `p.t_save` as a subset of the original time axis `t` to reduce the number of time steps that appear in the output (to avoid the output becoming large). If `p.t_save` is not specified, the model will output all time steps.

- `t` is the time axis for the simulation, in units of days.
    - `t(1)` is the start of the simulation.
    - `t(end)` is the end of the simulation.
    - The model time step is calculated as `t(i+1) - t(i)` (so doesn't have to be the same at every time step).
    
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

## Examples

### Example 1 - subglacial discharge

### Example 2 - intermediary circulation driven by shelf variability

### Example 3 - icebergs

### Example 4 - combined simulation
