function s = compute_fluxes(i, p, f, s)

% COMPUTE_FLUXES compute fluxes in the zmodel simulation.
%   S = COMPUTE_FLUXES(I, P, F, S) calls functions to compute the
%   plume fluxes, shelf fluxes, mixing fluxes, iceberg fluxes and vertical
%   fluxes for parameters structure P, boundary conditions F, and solution
%   S at timestep I.

% Calculate plume fluxes.
[s.QVg(:,i), s.QTg(:,i), s.QSg(:,i)] = get_plume_fluxes(i, p, f, s);

% Calculate shelf fluxes, return fluxes and exterior variables.
[s.QVs(:,i), s.QTs(:,i), s.QSs(:,i), s.Se(:,i), s.Te(:,i), s.phi(:,i)] = get_shelf_fluxes(i, p, f, s);

% Calculate tracer vertical mixing fluxes.
[s.QVk(:,i), s.QTk(:,i), s.QSk(:,i)] = get_mixing_fluxes(i, p, s, s.QVg(:,i), s.QVs(:,i));

% Calculate iceberg fluxes.
[s.QVi(:,i), s.QTi(:,i), s.QSi(:,i)] = get_iceberg_fluxes(i, p, s);

% Calculate vertical fluxes.
[s.QVv(:,i), s.QTv(:,i), s.QSv(:,i)] = get_vertical_fluxes(i, s, s.QVg(:,i)-s.QVs(:,i)+s.QVi(:,i));

end
