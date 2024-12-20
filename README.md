
# Fjord Reduced Physics Model (FjordRPM)

The Fjord Reduced Physics Model (FjordRPM) is an efficient "1.5-dimensional" model for simulating the dynamics of glacial fjords. The model splits a fjord into a number of vertically-stacked layers and solves for the evolution of layer properties due to processes that currently include (i) plumes driven by subglacial discharge, (ii) exchange with the continental shelf, (iii) iceberg melt and upwelling, (iv) vertical mixing and advection within the fjord. The model can represent multiple plumes and also the presence of a sill. Full details on the physics, numerical implementation and validation of the model can be found in the model description paper. In particular, the appendix of the model description paper contains an overview of all model parameters and variables.

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
    
- `a` is a structure containing the desired layer thicknesses, the initial conditions and the iceberg surface area profile. 
    - `a.H0` contains the desired layer thickness (m), which do not need to be uniform. They must sum to the depth of the fjord. The word "desired" is used here because the model will adjust the thicknesses slightly to ensure that the sill depth coincides with a layer boundary, making the treatment of the shelf exchange fluxes simpler. The actual layer thicknesses used are output in the solution structure as `s.H`.
    - `a.T0` contains the initial conditions for fjord temperature ($^{\circ}$C) - i.e., the initial temperature of the layers having thickness `a.H0`.
    - `a.S0` contains the initial conditions for fjord salinity (unitless).
    - `a.I0` contains the surface area of icebergs (m$^2$) in each of the layers having thickness `a.H0`.
    - All initial conditions should have dimensions `[p.N,1]`
    
- `s` is the output structure containing all solution variables. There are many such variables and we only summarise and pick out some key ones here, grouped by dimension.
    - The variable `s.t` contains the time axis for the solution (`s.t=p.t_save` if `p.t_save` is specified).
    - Variables of dimension `[p.N,1]` are layer variables that are static through the simulation, including the layer thickness `s.H`, layer iceberg surface area `s.I` and depth of the midpoint of each layer `s.z`, the last of which is useful for plotting.
    - Variables of dimension `[p.N,length(s.t)]` are time- and layer-dependent variables. These include the layer temperature `s.T` and salinity `s.S`, and the layer iceberg melt rate `s.mi`. The exchange fluxes are denoted `s.QXy` where X is the quantity being exchanged (V=volume, T=heat, S=salt) and y is the process doing the exchanging (s=shelf, i=icebergs, k=vertical mixing, v=vertical advection). Also included are the shelf temperature `s.Ts` and salinity `s.Ss`, which are the same as the input forcings but interpolated onto the model layers and output time steps. 
    - All plume-related solution variables have dimension `[n_plumes,p.N,length(s.t)]`, including for example the plume volume exchange fluxes `s.QVp` and plume submarine melt rate `s.mp`.
    - The subglacial discharge `s.Qsg` is the same as the input forcing `f.Qsg` but interpolated onto the output time steps.
    
**Note 1:** if a directory string is provided, for example `s = run_model(p,t,f,a,'directory_path/filename.mat')`, then all inputs and outputs will be saved in `directory_path/filename.mat`.

**Note 2:** full lists of model parameters and variables are provided in the appendix of the model description paper.

## Examples

In the `examples` directory we provide 4 idealised examples illustrating the various capabilities of the model and providing a basis to start from for building your own simulations. Each contains a *.m* file that runs the example and is heavily commented with explanation, and a *.mat* file that contains the result of the running the simulation. At the end, the scripts call `plotting_code/plotrpm` to make basic plots of the simulation. There is also a line that calls `plotting_code/animate` to create the mp4 animations that are included in the `examples` directory, but this line is commented as it takes a few minutes to create the animation.

### Example 1 - subglacial discharge

This is a simulation of the fjord response to seasonal variability in subglacial discharge. The fjord is 60 km long, 6 km wide, 800 m deep and has a sill of depth 400 m. There is a single glacier with grounding line depth 800 m and a single plume of width 250 m (though there is commented code to have multiple plumes/glaciers). The shelf conditions are uniform in depth and time with a temperature of 3 $^{\circ}$C and salinity of 34. Within each year, the subglacial discharge is Gaussian in time with a peak of 300 m$^3$/s on the 200th day of the year. There are no icebergs. The simulation runs for 3 years with a time step of 0.2 days and 40 model layers of equal, 20 m, thickness. The output is saved once a day.

The solution (run simulation to generate plots or see `examples/example1_subglacial_discharge.mp4`) shows the subglacial discharge-driven plume reaching the surface or close to the surface before intruding into the fjord. This intrusion results in seasonal cooling and freshening of the top 100 m of the fjord, which sets up an exchange with the shelf over the sill, whereby the top ~100 m flows from fjord to shelf, reaching a maximum of 0.15 m/s, and the ~100-400 m flows from shelf to fjord. When the subglacial discharge is active, there is some downwelling in the fjord where deep waters entrained by the plume have to be replaced by waters coming in over the sill. There is also seasonal plume-induced submarine melting reaching 6 m/d.

### Example 2 - intermediary circulation driven by shelf variability

### Example 3 - icebergs

### Example 4 - combined simulation
