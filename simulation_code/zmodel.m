function s = zmodel(p, t, f, a, path_out)

% ZMODEL z-model simulation.
%   [S, F] = ZMODEL(P, T, F, A, PATH_OUT) runs the z-model simulation for
%   parameters structure P, time T, forcings structure F, initial
%   conditions structure A, and returns solution structure S which includes
%   input forcing in the same time steps as S. If PATH_OUT is specified,
%   ZMODEL will save a file.

%% Check for errors in the given boundary and initial conditions.
[status, a] = check_inputs(p, a, t);

%% Preallocate and initialise variables- fluxes and box tracers-
% according to the number of layers for each timestep and store in s.
s = initialise_variables(p, a, t);

%% The main ZMODEL loop.
for i = 1:length(t)-1

    % Homogenise the heat and salt content of layers if the density
    % stratification is unstable and store the updated solution in s.
    s = homogenise_unstable_layers(i, p, s);

    % Compute the fluxes Q at the boundaries of each layer at timestep i
    % and associated exterior E and ice I variables.
    [Q, E] = compute_fluxes(i, p, f, s);
    
    % Step the fjord forwards to compute the tracer variables at timestep
    % i+1.
    Tr = step_solution_forwards(i, p, s, Q);

    % Store the fluxes computed at timestep i and the zmodel variables
    % computed at timestep i+1 in s.
    s = store_solution(i, s, Q, E, Tr);

    % Optional runtime plotting (for debugging).
    if p.plot_runtime
        plot_runtime_profile(i, t, f, p, s);
    end

    % Check for errors in the layer depths in the new solution.
    status = check_solution(p, s.H(:,i+1));
    
end

% Compute and store the fluxes at the final timestep.
s = compute_final_fluxes(p, f, s, t);

%% Get the output solution (save daily values to avoid large output files).
s = get_final_output(p, f, t, s, status);

%% Save output if a path and file name are provided.
if nargin > 4
    save_output(s, f, t, p, a, path_out);
end

end