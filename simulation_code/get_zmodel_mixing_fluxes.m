function [QVk0,QTk0,QSk0] = get_zmodel_mixing_fluxes(i, p, s)

% GET_ZMODEL_MIXING_FLUXES Compute mixing fluxes for the zmodel.
%   [QVK0, QTK0, QSK0] = GET_ZMODEL_MIXING_FLUXES(I, P, F, S) computes the
%   mixing fluxes QVK0, QTK0, QSSK0 for the given parameters P, boundary
%   conditions F and solution S at timestep I.

% Get tracer variables at timestep i.
H0 = s.H(:,i); T0 = s.T(:,i); S0 = s.S(:,i);

% The net volume mixing fluxes are always zero.
QVk0 = 0*H0;

% Compute the mixing fluxes going in/out of each layer; these are zero if
% the controlling parameter p.K0 is zero.
QS = 2*p.W*p.L*p.K0*(S0(2:end)-S0(1:end-1))./(H0(2:end)+H0(1:end-1));
QT = 2*p.W*p.L*p.K0*(T0(2:end)-T0(1:end-1))./(H0(2:end)+H0(1:end-1));

% The final fluxes are the sum of the fluxes going in and out of each
% layer.
QTk0 = [QT;0]-[0;QT];
QSk0 = [QS;0]-[0;QS];

end