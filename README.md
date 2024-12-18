
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
    - The model time step is calculated as `t(i+1) - t(i)` (so the time step can vary through the simulation if desired).
    
- `f` is a structure containing all the model forcings.
    - `f.ts` is a vector, dimensions `[1,length(f.ts)]`, of times (in days) on which the shelf properties are defined.
    - `f.zs` is a vector, dimensions `[length(f.zs),1]`, of depths (m, negative below sea level) on which the shelf properties are defined.
    - `f.Ts` is an array, dimensions `[length(f.zs),length(f.ts)]`, of the shelf temperature ($^{\circ}$C) at the times `f.ts` and depths `f.zs`.
    - `f.Ss` is an array, dimensions `[length(f.zs),length(f.ts)]`, of the shelf salinity (unitless) at the times `f.ts` and depths `f.zs`.
    - `f.tsg` is a vector, dimensions `[1,length(f.tsg)]`, of times (in days) at which the subglacial discharge is defined.
    - `f.Qsg` is an array, dimensions `[n_plumes,length(f.tsg)]`, of the subglacial discharge (m$^3$/s) at the times `f.tsg` for each of the number `n_plumes` of subglacial discharge-driven plumes.
    - If `n_plumes > 1`, then `p.Hgl` and `p.Wp` both need to be arrays containing the potentially different grounding line depths and plume widths for each plume (for an example, see comments in `examples/example1_subglacial_discharge`).
    - `f.ts` and `f.tsg` do not need to be the same as each other or the same as `t`, but they do need to all be consistent and the forcing needs to specified for the full duration of the simulation. Similarly `f.zs` does not need to coincide with the model layers. When the time/depth axes for the specification of the forcing do not coincide with the model time stepping and layers, the forcings are interpolated in time and depth to coincide with the model time stepping and layers.
    
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
