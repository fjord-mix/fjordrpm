function [s,f] = zmodel(p,t,f,a,path_out)

% ZMODEL z-model simulation.
%   [S,F] = ZMODEL(P, T, F, A, PATH_OUT) runs the z-model simulation for parameters structure P,
%   time T, forcings structure F, initial conditions structure A, and returns solution structure
%   S and forcing structure F in the same time steps as S. If PATH_OUT is
%   specified, ZMODEL will save a file. If F and/or A are not specified, default idealised values will be used.

%% Check boundary and initial conditions.
% Set idealised boundary and initial conditions, if not given, based on input parameters.
% Boundary conditions:
if nargin < 3, f = get_idealised_forcing(p, t); end
if isempty(f), f = get_idealised_forcing(p, t); end % we cannot use an OR statement here
% Initial conditions:
if nargin < 4, a = get_initial_conditions(p, f); end
if isempty(a), a = get_initial_conditions(p, f); end

% Check for errors in the given initial state.
status = check_zmodel_initialisation(p, a);

% Preallocate and initialise variables- fluxes and box tracers-
% according to the number of layers for each timestep and store in s.
s = initialise_zmodel_variables(p, f, a, t);

%% Optional runtime plotting (for debugging)
if p.plot_runtime
    % hf_track = monitor_boxmodel([],1,H,T,S,f);
    % hf_track = show_boxmodel([],1,t,H,T,S,[],[],[],[],f);
    % s_bnds = [min(f.Ss(:)) max(f.Ss(:))+0.1];
    % plot_debug_profile(0,t,f,p,a.H0',S,s_bnds,T,[]);
end

%% The main ZMODEL loop.
for i = 1:length(t)-1

    % Homogenise the heat and salt content of layers if the density
    % stratification is unstable and store the updated solution in s.
    s = homogenise_zmodel_unstable_layers(i, p, s);

    % Compute the fluxes at the boundaries of each layer at timestep i.
    Q_i = compute_zmodel_fluxes(i, p, f, s);
    
    % Step the fjord forwards to compute the zmodel variables at timestep i+1.
    Z_ip1 = step_zmodel_forwards(i, p, s, Q_i);

    % Store the fluxes computed at timestep i and the zmodel variables computed at
    % timestep i+1 in s.
    s = store_zmodel_solution(i, s, Q_i, Z_ip1);

    % Optional runtime plotting (for debugging).
    if p.plot_runtime
        % hf_track = monitor_boxmodel(hf_track,i,H,T,S,f);
        % hf_track = show_boxmodel([],i,t,H,T,S,QVs,QVg,QVk,QVb,f);
        % plot_debug_profile(i,t,f,p,H,S,[],T,[]);
        plot_runtime_profile(i,t,f,p,s.H,s.S,s.T,s.Se,s.QVs)
    end

    % Check for errors in the new state.
    status = check_zmodel(p, s.H(:,i+1));
    
end

%% Get the output solution (save daily values to avoid large output files).
s = get_zmodel_output(p, f, t, s, status);

%% Save output if a path+file name are provided
if nargin > 4
    fjord_output.s = s;
    fjord_output.f = f;
    fjord_output.t = t;
    fjord_output.p = p;
    fjord_output.a = a;
    save(path_out,'fjord_output','-v7.3'); % v7.3 allows large files (> 2GB), which might happen in very long runs
end

end