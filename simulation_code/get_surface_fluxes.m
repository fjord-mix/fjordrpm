function [QVsurf0, QTsurf0, QSsurf0] = get_surface_fluxes(i, p, s)

% GET_SURFACE_FLUXES Compute surface fluxes.
%   [QVa0, QTa0, QSa0] = GET_SURFACE_FLUXES(i, p, s)
%   computes the surface fluxes QVa0, QTa0, QSa0
%   for the given parameters p at timestep i.

% get tracer variables at timestep i
H0 = s.H; T0 = s.T(:,i); S0 = s.S(:,i);

% initialise outputs
[QVsurf0, QTsurf0, QSsurf0] = deal(0*H0);

% riverine input at timestep i goes into top layer
QVsurf0(1) = s.Qr(:,i);

% heat and salt fluxes remain 0 because T=0, S=0 assumed for riverine input

end