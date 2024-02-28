function [QIi0,QTi0,QSi0,M0,QVmi0,QTmi0,QSmi0] = get_iceberg_fluxes(H0,T0,S0,I0,zi,p)


% GET_ICBERG_FLUXES Compute iceberg fluxes.
%   [QII0,QTI0,QSI0,M0] = GET_ICEBERG_FLUXES(H0,T0,S0,I0,ZI,P)
%   computes the iceberg fluxes for the given parameters.

if p.M0==0 % if no icebergs
    QIi0 = 0*H0;
    QTi0 = 0*H0;
    QSi0 = 0*H0;
    QVmi0 = 0*H0;
    QTmi0 = 0*H0;
    QSmi0 = 0*H0;
    M0 = 0*zi;

else

    % melting in boxes
    ints = [0;-cumsum(H0)];
    zj = 0.5*(ints(1:end-1)+ints(2:end)); % mean depth of box
    Tf = p.l1*S0 + p.l2 + p.l3*zj; % local freezing point

    % get vector of iceberg concentration that resolves box boundaries
    zi0 = unique(sort([0,zi,-cumsum(H0)']));
    Ii0 = interp1(zi,I0,zi0,'pchip','extrap');

    % Preallocate variables
    QIi0 = zeros(1, length(ints)-1);
    % do numerical integral to get melt fluxes
    for k=1:length(ints)-1
        inds = find(zi0<=ints(k) & zi0>=ints(k+1));
        QIi0(k) = p.W*p.L*p.M0*(T0(k)-Tf(k))*trapz(zi0(inds),Ii0(inds));
    end
    QTi0 = -QIi0*p.l/p.cw;
    QSi0 = -QIi0.*S0';

    % melt profile defined on zi reference depths
    M0 = 0*zi;
    for k=1:length(ints)-1
        inds = zi<=ints(k) & zi>=ints(k+1);
        M0(inds) = p.M0*(T0(k)-Tf(k));
    end

    % mixing fluxes between boxes to account for buoyant rising of
    % freshwater
    % Preallocate variables
    [QVmI, QTmI, QSmI, gmelt, Heff, SA_ice] =  deal(zeros(1, p.N+p.sill-1));
    % loop over interfaces, allow plumes to rise from below sill layer
    for k=1:p.N+p.sill-1
        % reduced gravity
        gmelt(k) = p.g*(p.betaS*S0(k)-p.betaT*T0(k)); % relative to the iceberg
        % effective height rise of plume in layer
        inds = find(zi0<=ints(k) & zi0>=ints(k+1));
        Heff(k) = trapz(zi0(inds),Ii0(inds));
        % Surface area of icebergs to melt (assume characteristic length scale 400)
        SA_ice(k) = p.W*p.L*6*sqrt(6)/(p.H/2)*Heff(k);
        % fluxes
        if Heff(k) == 0 % no iceberg concentration (e.g. on first timestep)
            QVmI(k) = 0;
            QTmI(k) = 0;
            QSmI(k) = 0;
        elseif gmelt(k) < 0 % iceberg plume is denser than the box
            QVmI(k) = 0;
            QTmI(k) = 0;
            QSmI(k) = 0;
        else % iceberg plume which mixes into box above 
            QVmI(k) = p.gamma*QIi0(k+1)*p.alphaI^(2/3)*gmelt(k)^(1/3)*(SA_ice(k)/Heff(k))^(2/3)*Heff(k)/2;
            QTmI(k) = QVmI(k)*T0(k+1)+p.gamma*QTi0(k+1);
            QSmI(k) = QVmI(k)*S0(k+1)+p.gamma*QSi0(k+1);
        end

    end

    % final fluxes
    QVmi0 = [QVmI,0]'-[0,QVmI]';
    QTmi0 = [QTmI,0]'-[0,QTmI]';
    QSmi0 = [QSmI,0]'-[0,QSmI]';

    % if p.sill % ensures there is no mixing at the sill layer interface
    %     QVmi0(p.N+p.sill) = 0;
    %     QTmi0(p.N+p.sill) = 0;
    %     QSmi0(p.N+p.sill) = 0;
    % end


end

end