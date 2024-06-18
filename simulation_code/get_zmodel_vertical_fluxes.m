function [QVv0,QTv0,QSv0] = get_zmodel_vertical_fluxes(i, p, s, Qnet)

% GET_ZMODEL_VERTICAL_FLUXES Compute vertical fluxes for the zmodel.
%   [QVV0, QTV0, QSV0] = GET_ZMODEL_VERTICAL_FLUXES(I, P, S, QNET) computes
%   the shelf fluxes QVV0, QTV0, QSV0 for the solution S at timestep I with
%   net flux QNET.

% Get tracer variables at timestep i.
T0 = s.T(:,i); S0 = s.S(:,i);

% The vertical flux required for no net volume change is the sum of the flux imbalances above
QVint = -cumsum(Qnet(1:end-1));
% The relevant temperature/salinity flux depends on the direction.
QTint = (QVint<0).*QVint.*T0(1:end-1) + (QVint>=0).*QVint.*T0(2:end);
QSint = (QVint<0).*QVint.*S0(1:end-1) + (QVint>=0).*QVint.*S0(2:end);

% The final fluxes are the sum of the fluxes going in and out of each
% layer.
QVv0 = [QVint;0]-[0;QVint];
QTv0 = [QTint;0]-[0;QTint];
QSv0 = [QSint;0]-[0;QSint];

end