function [QVsurf0, QTsurf0, QSsurf0] = get_surface_fluxes(i, p, s)

% GET_SURFACE_FLUXES Compute surface fluxes.
%   [QVa0, QTa0, QSa0] = GET_SURFACE_FLUXES(i, p, s)
%   computes the surface fluxes QVa0, QTa0, QSa0
%   for the given parameters p at timestep i.

% get tracer variables at timestep i
H0 = s.H; T0 = s.T(:,i); S0 = s.S(:,i);

% initialise outputs
[QVsurf0, QTsurf0, QSsurf0] = deal(0*H0);

% surface forcings at timestep i
Qr0 = s.Qr(i); % river volume
Tr0 = s.Tr(i); % river temperature
Sr0 = s.Sr(i); % river salinity
Ta0 = s.Ta(i); % air temperature

% fluxes - riverine and air-sea heat flux
QVsurf0(1) = Qr0;
QTsurf0(1) = Qr0*Tr0 + (p.kairsea/(p.rhoref*p.cw))*p.W*p.L*(s.Ta(i)-T0(1));
QSsurf0(1) = Qr0*Sr0;

end