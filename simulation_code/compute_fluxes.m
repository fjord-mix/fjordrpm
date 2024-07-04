function [Q, E] = compute_fluxes(i, p, f, s)

% COMPUTE_FLUXES compute fluxes in the zmodel simulation.
%   [Q, E, I] = COMPUTE_FLUXES(I, P, F, S) calls functions to compute the
%   plume fluxes, shelf fluxes, mixing fluxes, iceberg fluxes and vertical
%   fluxes for parameters structure P, boundary conditions F, and solution
%   S at timestep I and returns fluxes Q, exterior variables E and iceberg
%   variables I.

% Calculate plume fluxes.
[Q.QVg, Q.QTg, Q.QSg] = get_plume_fluxes(i, p, f, s);

% Calculate shelf fluxes, return fluxes and exterior variables.
[Q.QVs, Q.QTs, Q.QSs, E.Se, E.Te, E.phi] = get_shelf_fluxes(i, p, f, s);

% Calculate tracer vertical mixing fluxes.
[Q.QVk, Q.QTk, Q.QSk] = get_mixing_fluxes(i, p, s, Q.QVg, Q.QVs);

% Calculate iceberg fluxes.
[Q.QVi, Q.QTi, Q.QSi] = get_iceberg_fluxes(i, p, s);

% Calculate vertical fluxes.
[Q.QVv, Q.QTv, Q.QSv] = get_vertical_fluxes(i, s, Q.QVg-Q.QVs+Q.QVi);

end
