function [QVk0,QTk0,QSk0] = get_mixing_fluxes(H0,T0,S0,p)

% GET_MIXING_FLUXES Compute vertical mixing fluxes.
%   [QVk0,QTk0,QSk0] = GET_MIXING_FLUXES(H0,T0,S0,P)
%   computes the vertical mixing fluxes for the given parameters.

% net volume fluxes are always 0
QVk0 = 0*H0;

% if no vertical mixing
if p.K0==0  
    QTk0 = 0*H0;
    QSk0 = 0*H0;

else
    % Preallocate variables
    [Q, QT, QS] =  deal(zeros(1, p.N+p.sill-1));

    % loop over interfaces
    for k=1:p.N+p.sill-1
        QS(k) = 2*p.W*p.L*p.K0*(S0(k+1)-S0(k))/(H0(k+1)+H0(k));
        QT(k) = 2*p.W*p.L*p.K0*(T0(k+1)-T0(k))/(H0(k+1)+H0(k));
    end

    % final fluxes
    QTk0 = [QT,0]'-[0,QT]';
    QSk0 = [QS,0]'-[0,QS]';

end

end