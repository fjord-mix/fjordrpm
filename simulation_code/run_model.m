function s = run_model(p, t, f, a, path_out)

% RUN_MODEL Run fjordRPM simulation.
%   s = RUN_MODEL(p, t, f, a, path_out) runs the fjordRPM simulation for
%   parameters structure p, time vector t, forcings structure f, initial
%   conditions structure a, and returns solution structure s. If path_out
%   is specified, RUN_MODEL will save a file.

%% Check for errors in the given inputs
status = check_inputs(p, t, f, a);

%% Preallocate and initialise variables
s = initialise_variables(p, t, f, a);

%% The timestepping loop
for i = 1:length(t)-1

    % Homogenise heat & salt content of layers where density is unstable
    s = homogenise_unstable_layers(i, p, s);
    
    % Compute the fluxes ready for next timestep
    s = compute_fluxes(i, p, s);

    % Step the tracers forward in time
    s = step_solution_forwards(i, p, s);

    % Optional runtime plotting
    if p.plot_runtime
        plot_runtime_profile(i, p, t, s);
    end

    % Check for errors at this new timestep
    status = check_solution(i, s);
    
end

%% Subsample solution structure to requested output frequency
s = get_final_output(p, t, s, status);

%% Save output if a path and file name are provided
if nargin > 4
    save_output(p, t, f, a, s, path_out);
end

end