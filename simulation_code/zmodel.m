function [s,f] = zmodel(p,t,f,a,path_out)

% ZMODEL z-model simulation.
%   [S,F] = ZMODEL(P,T,F,A,PATH_OUT) runs the z-model simulation for parameters structure P,
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

% Preallocate and initialise variables- fluxes and box tracers- according to the number of layers for each timestep.
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
    % stratification is unstable.
    [s.T(:,i), s.S(:,i), s.VT(:,i), s.VS(:,i)] = ...
        homogenise_unstable_layers(p, s.V(:,i), s.T(:,i), s.S(:,i), s.VT(:,i), s.VS(:,i));

     % Compute the fluxes at the boundaries of each layer.
    [s.QVg(:,i),s.QTg(:,i),s.QSg(:,i),...
        s.QVs(:,i),s.QTs(:,i),s.QSs(:,i),s.Se(:,i),s.Te(:,i),s.phi(:,i),...
        s.QVk(:,i),s.QTk(:,i),s.QSk(:,i),...
        s.QIi(:,i),s.QTi(:,i),s.QSi(:,i),s.M(:,i),s.QVmi(:,i),s.QTmi(:,i),s.QSmi(:,i),...
        s.QVv(:,i),s.QTv(:,i),s.QSv(:,i)] ...
        = compute_fluxes(...
        s.H(:,i),s.T(:,i),s.S(:,i),f.Qsg(i),p,f.zs,f.Ts(:,i),f.Ss(:,i), ...
        s.V(:,i),s.I(:,i),f.zi, s);
    
    % Step the fjord forwards.
    [s.V(:,i+1), s.VT(:,i+1), s.VS(:,i+1), s.H(:,i+1), s.T(:,i+1), s.S(:,i+1), s.I(:,i+1)] = ...
        step_zmodel_forwards(p, f, s.V(:,i), s.VT(:,i), s.VS(:,i), s.I(:,i), s.M(:,i), f.D(i),...
        s.QVg(:,i), s.QVs(:,i), s.QVk(:,i), s.QVmi(:,i), s.QVv(:,i), ...
        s.QTg(:,i), s.QTs(:,i), s.QTk(:,i), s.QTi(:,i), s.QTmi(:,i), s.QTv(:,i), ...
        s.QSg(:,i), s.QSs(:,i), s.QSk(:,i), s.QSi(:,i), s.QSmi(:,i), s.QSv(:,i));

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
s = get_zmodel_output(p, f, t, status, s.H, s.T, s.S, s.V, s.I, ...
    s.QVg, s.QTg, s.QSg, s.QVs, s.QTs, s.QSs, s.Se, s.Te, s.phi, s.QVk, s.QTk, s.QSk, ...
    s.QVv, s.QTv, s.QSv, s.QIi, s.QTi, s.QSi, s.QVmi, s.M);


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