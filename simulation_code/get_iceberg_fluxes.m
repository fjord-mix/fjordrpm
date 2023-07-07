function [QIi0,QTi0,QSi0,M0] = get_iceberg_fluxes(H0,T0,S0,I0,zi,p)

% GET_ICBERG_FLUXES Compute iceberg fluxes.
%   [QII0,QTI0,QSI0,M0] = GET_ICEBERG_FLUXES(H0,T0,S0,I0,ZI,P)
%   computes the iceberg fluxes for the given parameters.

if p.M0==0 % if no icebergs
    QIi0 = 0*H0;
    QTi0 = 0*H0;
    QSi0 = 0*H0;
    M0 = 0*zi;

else
    % melting in boxes
    ints = [0;cumsum(H0)];
    zj = 0.5*(ints(1:end-1)+ints(2:end)); % mean depth of box
    Tf = p.l1*S0 + p.l2 + p.l3*zj; % local freezing point

    % get vector of iceberg concentration that resolves box boundaries 
    zi0 = unique(sort([zi,-cumsum(H0)']));
    Ii0 = interp1(zi,I0,zi0,'pchip','extrap');

    % Preallocate variables
    QIi0 = zeros(1, length(ints)-1);
    % do numerical integral to get melt fluxes
    for k=1:length(ints)-1
        inds = find(zi0<=ints(k) & zi0>=ints(k+1));
        QIi0(k) = p.W*p.L*p.M0*(T0(k)-Tf(k))*trapz(zi0(inds),Ii0(inds)); % modified Tf0 -> Tf
    end
    QTi0 = -QIi0*p.l/p.cw;
    QSi0 = -QIi0.*S0'; % modified S0 -> S0'

    % melt profile defined on zi reference depths
    M0 = 0*zi;
    for k=1:length(ints)-1
        inds = zi<=ints(k) & zi>=ints(k+1);
        M0(inds) = p.M0*(T0(k)-Tf(k));
    end
end

end