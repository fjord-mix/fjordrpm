function s = zmodel(p, t, f, a, path_out)

% ZMODEL z-model simulation.
%   [S, F] = ZMODEL(P, T, F, A, PATH_OUT) runs the z-model simulation for
%   parameters structure P, time T, forcings structure F, initial
%   conditions structure A, and returns solution structure S which includes
%   input forcing in the same time steps as S. If PATH_OUT is specified,
%   ZMODEL will save a file.

%% Check for errors in the given boundary and initial conditions
[status, a] = check_inputs(p, a, t);

%% Preallocate and initialise variables
s = initialise_variables(p, t, f, a);

%% The timestepping loop
for i = 1:length(t)-1

    % Homogenise heat & salt content of layers where density is unstable
    s = homogenise_unstable_layers(i, p, s);
    
    % Compute the fluxes ready for next timestep
    s = compute_fluxes(i, p, f, s);

    % Step the tracers forward in time
    s = step_solution_forwards(i, p, s);

    % Optional runtime plotting
    if p.plot_runtime
        plot_runtime_profile(i, t, f, p, s);
    end

    % Check for errors at this new timestep
    status = check_solution(i, s);
    
end

%% Subsample solution structure to requested output frequency
s = get_final_output(p, f, t, s, status);

%% Save output if a path and file name are provided
if nargin > 4
    save_output(s, f, t, p, a, path_out);
end

end