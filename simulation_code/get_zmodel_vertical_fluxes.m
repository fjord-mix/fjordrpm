function Qv = get_zmodel_vertical_fluxes(i, p, s, Q)

% GET_VERTICAL_FLUXES Compute vertical fluxes when thicknesses fixed.
%   [QVv0,QTv0,QSv0] = GET_VERTICAL_FLUXES(Qnet,T0,S0,p)
%   computes the vertical fluxes for the given parameters.

Qnet = Q.Qg.V-Q.Qs.V+Q.Qi.Vm;
T0 = s.T(:,i);
S0 = s.S(:,i);

% if fixed thickness layers
if p.fixedthickness

    for i=1:length(Qnet)-1
        % vertical flux required for no net volume change is the sum
        % of the flux imbalances above
        QVint(i) = -sum(Qnet(1:i));
        % relevant T/S for flux depends on direction
        if QVint(i)<0
            QTint(i) = QVint(i)*T0(i);
            QSint(i) = QVint(i)*S0(i);
        else
            QTint(i) = QVint(i)*T0(i+1);
            QSint(i) = QVint(i)*S0(i+1);
        end
    end
    
    % final fluxes
    QVv0 = [QVint,0]'-[0,QVint]';
    QTv0 = [QTint,0]'-[0,QTint]';
    QSv0 = [QSint,0]'-[0,QSint]';

% if variable thickness layers these fluxes are all 0
else

    QVv0 = zeros(1,p.N);
    QTv0 = zeros(1,p.N);
    QSv0 = zeros(1,p.N);

end
Qv.V = QVv0;
Qv.T = QTv0;
Qv.S = QSv0;
end