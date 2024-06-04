function [s,f] = zmodel(p,t,f,a,path_out)

% ZMODEL One-dimensional z-model simulation.
%   [S,F] = ZMODEL(P,T,F,A,PATH_OUT) runs the z-model simulation for parameters structure P,
%   time T, forcings structure F, initial conditions structure A. Returns solution structure
%   S and forcing structure F in the same time steps as S. If PATH_OUT is
%   specified, will save a file.
%   If F and/or A are not specified, default idealised values will be used.

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

%% Preallocate and initialise variables according to the number of layers for each timestep.
[V, T, S, I, VT, VS, ...
    QVg, QTg, QSg,...
    QVs, QTs, QSs,...
    QVk, QTk, QSk,...
    QVmi, QTmi, QSmi,...
    QIi, QTi, QSi,...
    QVv, QTv, QSv,...
    Te, Se, phi, M] = initialise_zmodel_variables(p, f, a, t);

%% Optional runtime plotting (for debugging)
if p.plot_runtime
    % hf_track = monitor_boxmodel([],1,H,T,S,f);
    % hf_track = show_boxmodel([],1,t,H,T,S,[],[],[],[],f);
    plot_runtime_profile(0,t,f,p,a.H0,S,T,[]);
end

%% The main ZMODEL loop.
for i = 1:length(t)-1

    % Homogenise the heat and salt content of layers if the density
    % stratification is unstable.
    [T(:,i), S(:,i), VT(:,i), VS(:,i)] = ...
        homogenise_unstable_layers(p, V, T(:,i), S(:,i), VT(:,i), VS(:,i));

    % Compute the fluxes at the boundaries of each layer.
    [QVg(:,i),QTg(:,i),QSg(:,i),...
        QVs(:,i),QTs(:,i),QSs(:,i),Se(:,i),Te(:,i),phi(:,i),...
        QVk(:,i),QTk(:,i),QSk(:,i),...
        QIi(:,i),QTi(:,i),QSi(:,i),M(:,i),QVmi(:,i),QTmi(:,i),QSmi(:,i),...
        QVv(:,i),QTv(:,i),QSv(:,i)] ...
        = compute_fluxes(...
        a.H0', T(:,i),S(:,i),f.Qsg(i),p,f.zs,f.Ts(:,i),f.Ss(:,i), ...
        V,I(:,i),f.zi);

    % Step the fjord forwards.
    [VT(:,i+1), VS(:,i+1), T(:,i+1), S(:,i+1), I(:,i+1)] = ...
        step_zmodel_forwards(p, f, V, VT(:,i), VS(:,i), I(:,i), M(:,i), f.D(i),...
        QTg(:,i), QTs(:,i), QTk(:,i), QTi(:,i), QTmi(:,i), QTv(:,i), ...
        QSg(:,i), QSs(:,i), QSk(:,i), QSi(:,i), QSmi(:,i), QSv(:,i));

    % Plot model evolution (for optional runtime debugging).
    if p.plot_runtime
        % hf_track = monitor_boxmodel(hf_track,i,H,T,S,f);
        % hf_track = show_boxmodel([],i,t,H,T,S,QVs,QVg,QVk,QVb,f);
        plot_runtime_profile(i,t,f,p,H,S,T,[]);
    end

end

%% Get the output solution (save daily values to avoid large output files).
s = get_zmodel_output(p, f, t, status, T, S, I, ...
    QVg, QTg, QSg, QVs, QTs, QSs, Se, Te, phi, QVk, QTk, QSk, ...
    QVv, QTv, QSv, QIi, QTi, QSi, QVmi, M);


%% Save output if a path+file name are provided.
if nargin > 4
    fjord_output.s = s;
    fjord_output.f = f;
    fjord_output.t = t;
    fjord_output.p = p;
    fjord_output.a = a;
    save(path_out,'fjord_output','-v7.3'); % v7.3 allows large files (> 2GB), which might happen in very long runs
end

end
