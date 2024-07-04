function [QVs0,QTs0,QSs0,Se0,Te0,phi0] = get_shelf_fluxes(i, p, f, s)

% GET_SHELF_FLUXES Compute shelf fluxes for the zmodel.
%   [QVS0,QTS0,QSS0,SE0,TE0,PHI0] = GET_SHELF_FLUXES(I, P, F, S) computes
%   the shelf fluxes QVS0, QTS0, QSSS0 and shelf variables SE0, TE0, PHI0
%   for the given parameters P, boundary conditions F and solution S at
%   timestep I.

% Get tracer variables at timestep i.
H0 = s.H; T0 = s.T(:,i); S0 = s.S(:,i);
% Get boundary conditions at timestep i.
Qsg0 = f.Qsg(i); zs = f.zs; Ts = f.Ts(:,i); Ss = f.Ss(:,i);
% Initialise variables.
[phi0, QVs0] = deal(zeros(p.N, 1));

% Calculate the mean shelf T/S over box model layers above sill.
[Te0, Se0] = bin_ocean_profiles(Ts, Ss, zs, H0);

% Compute the fjord-to-shelf reduced gravity.
gp = p.g*(p.betaS*(S0(1:s.ksill)-Se0(1:s.ksill))-p.betaT*(T0(1:s.ksill)-Te0(1:s.ksill)));

% Calculate the potentials over above-sill layers.
phi0(1) = gp(1)*H0(1)/2;
for k=2:s.ksill
    phi0(k) = phi0(k-1)+gp(k-1)*H0(k-1)/2+gp(k)*H0(k)/2;
end

% Compute the above-sill fluxes before barotropic compensation; these are
% zero if the controlling parameter p.C0 = 0.
Q = p.C0*p.W*H0(1:s.ksill).*phi0(1:s.ksill)/p.L;

% Compute the above-sill fluxes after barotropic compensation, ensuring
% depth mean = Qsg0 when plume is turned on, and 0 when plume is turned
% off.
if p.P0==0
    Qsg0 = 0;
end
QVs0(1:s.ksill) = Q + H0(1:s.ksill)*(Qsg0-sum(Q))/sum(H0(1:s.ksill));

% Compute the resulting heat/salt fluxes.
QTs0 = (QVs0>=0).*QVs0.*T0 + (QVs0<0).*QVs0.*Te0;
QSs0 = (QVs0>=0).*QVs0.*S0 + (QVs0<0).*QVs0.*Se0;

end
