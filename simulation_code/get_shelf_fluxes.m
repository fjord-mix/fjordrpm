function [QVs0,QTs0,QSs0,phi0] = get_shelf_fluxes(i, p, s)

% GET_SHELF_FLUXES Compute shelf fluxes for the zmodel.
%   [QVs0,QTs0,QSs0,phi0] = GET_SHELF_FLUXES(i, p, s) computes
%   the shelf fluxes QVs0, QTs0, QSs0 and working variable phi0
%   for the given parameters p and solution s at
%   timestep i.

% get tracer variables at timestep i
H0 = s.H; T0 = s.T(:,i); S0 = s.S(:,i); Ts0 = s.Ts(:,i); Ss0 = s.Ss(:,i);
% initialise variables
[phi0, QVs0] = deal(zeros(p.N, 1));

% compute the fjord-to-shelf reduced gravity
gp = p.g*(p.betaS*(Ss0(1:s.ksill)-S0(1:s.ksill))-p.betaT*(Ts0(1:s.ksill)-T0(1:s.ksill)));

% calculate the potentials over above-sill layers
phi0(1) = gp(1)*H0(1)/2;
for k=2:s.ksill
    phi0(k) = phi0(k-1)+gp(k-1)*H0(k-1)/2+gp(k)*H0(k)/2;
end

% compute the above-sill fluxes before barotropic compensation; these are
% zero if the controlling parameter p.C0 = 0
Q = p.C0*p.W*H0(1:s.ksill).*phi0(1:s.ksill)/p.L;

% compute the above-sill fluxes after barotropic compensation, ensuring
% depth mean = Qsg0+submarine melt flux+riverine flux
QVs0(1:s.ksill) = Q - H0(1:s.ksill)*(s.Qr(i)+sum(s.Qsg(:,i))+sum(sum(s.QMp(:,:,i)))+sum(Q))/sum(H0(1:s.ksill));

% compute the resulting heat/salt fluxes
QTs0 = (QVs0>=0).*QVs0.*Ts0 + (QVs0<0).*QVs0.*T0;
QSs0 = (QVs0>=0).*QVs0.*Ss0 + (QVs0<0).*QVs0.*S0;

end
