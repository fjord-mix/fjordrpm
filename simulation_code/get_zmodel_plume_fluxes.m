function [QpV0 ,QpT0 ,QpS0 ] = get_zmodel_plume_fluxes(i, p, f, s)

% GET_PLUME_FLUXES Compute plume fluxes.
%   [QPV0,QPT0,QPS0] = GET_PLUME_FLUXES(H0,T0,S0,QSG0,P)
%   computes the plume fluxes for the given parameters.

H0 = s.H(:,i);
T0 = s.T(:,i);
S0 = s.S(:,i);
Qsg0 = f.Qsg(i);

if Qsg0==0 || p.P0==0 % i.e. if no plume
    QpV0 = 0*H0;
    QpT0 = 0*H0;
    QpS0 = 0*H0;
else
    % find box model layer containing grounding line
    ints = cumsum(H0);
    kgl = find(ints>=abs(p.zgl)-1e-6, 1 );

    % plume properties at grounding line
    k = kgl;
    Qp(k) = Qsg0;
    Sp(k) = 0;
    Tp(k) = 0;
    gp(k) = p.g*(p.betaS*(S0(k)-Sp(k))-p.betaT*(T0(k)-Tp(k)));

    % properties at first interface above grounding line
    % need special treatment because box might be partial if
    % grounding line does not coincide with box boundaries
    if kgl>1 % check there is a box above grounding line
        k = k-1;
        Qp(k) = Qp(k+1) + p.P0^(2/3)*Qp(k+1)^(1/3)*gp(k+1)^(1/3)*(abs(p.zgl)-ints(k));
        Tp(k) = (Qp(k+1)*Tp(k+1)+(Qp(k)-Qp(k+1))*T0(k+1))/Qp(k);
        Sp(k) = (Qp(k+1)*Sp(k+1)+(Qp(k)-Qp(k+1))*S0(k+1))/Qp(k);
        gp(k) = p.g*(p.betaS*(S0(k)-Sp(k))-p.betaT*(T0(k)-Tp(k)));
    end

    % apply to successive interfaces higher provided plume is still rising
    while gp(k)>0 && k>1
        k = k-1;
        Qp(k) = Qp(k+1) + p.P0^(2/3)*Qp(k+1)^(1/3)*gp(k+1)^(1/3)*H0(k+1);
        Tp(k) = (Qp(k+1)*Tp(k+1)+(Qp(k)-Qp(k+1))*T0(k+1))/Qp(k);
        Sp(k) = (Qp(k+1)*Sp(k+1)+(Qp(k)-Qp(k+1))*S0(k+1))/Qp(k);
        gp(k) = p.g*(p.betaS*(S0(k)-Sp(k))-p.betaT*(T0(k)-Tp(k)));
    end

    % now calculate the resulting box fluxes
    knb = find(gp<0); 
    if isempty(knb)
        knb=1; 
    end % find neutral buoyancy box
    % flux in boxes below grounding line and above neutral buoyancy are 0    
    inds = find((1:length(H0))>kgl | (1:length(H0))<knb);

      % initialise vectors
    [QpV0, QpT0, QpS0] = deal(zeros(p.N, 1));

    QpV0(inds) = 0;
    QpT0(inds) = 0;
    QpS0(inds) = 0;
    % flux in boxes from grounding line to below neutral buoyancy
    for k=kgl:-1:knb+1
        QpV0(k) = Qp(k)-Qp(k-1);
        QpT0(k) = QpV0(k)*T0(k);
        QpS0(k) = QpV0(k)*S0(k);
    end
    % flux into neutral buoyancy box
    QpV0(knb) = Qp(knb);
    QpT0(knb) = Qp(knb)*Tp(knb);
    QpS0(knb) = Qp(knb)*Sp(knb);
end

% Transpose vectors before outputting
% QpV0 = QpV0'; QpT0 = QpT0'; QpS0 = QpS0';

% Qg.V = QpV0;
% Qg.T = QpT0;
% Qg.S = QpS0;
% 
% s.QVg(:,i) = Q.Qg.V;
% s.QTg(:,i) = Q.Qg.T;
% s.QSg(:,i) = Q.Qg.S;

end
