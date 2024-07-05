function s = compute_fluxes(i, p, f, s)

% COMPUTE_FLUXES compute fluxes in the simulation.
%   s = COMPUTE_FLUXES(i, p, f, s) calls functions to compute the
%   plume fluxes, exterior fluxes, mixing fluxes, iceberg fluxes and vertical
%   fluxes for parameters structure p, boundary conditions f, and solution
%   s at timestep i.

% Calculate plume fluxes
[s.QVp(:,i), s.QTp(:,i), s.QSp(:,i)] = get_plume_fluxes(i, p, f, s);

% Calculate shelf fluxes
[s.QVs(:,i), s.QTs(:,i), s.QSs(:,i), s.Se(:,i), s.Te(:,i), s.phi(:,i)] = get_shelf_fluxes(i, p, f, s);

% Calculate tracer vertical mixing fluxes
[s.QVk(:,i), s.QTk(:,i), s.QSk(:,i)] = get_mixing_fluxes(i, p, s);

% Calculate iceberg fluxes
[s.QVi(:,i), s.QTi(:,i), s.QSi(:,i)] = get_iceberg_fluxes(i, p, s);

% Calculate vertical fluxes
[s.QVv(:,i), s.QTv(:,i), s.QSv(:,i)] = get_vertical_fluxes(i, s);

end
