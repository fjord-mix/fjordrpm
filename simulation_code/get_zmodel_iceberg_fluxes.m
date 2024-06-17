function [QIi0,QTi0,QSi0,M0,QVmi0,QTmi0,QSmi0] = get_zmodel_iceberg_fluxes(i, p, s)


% GET_ICEBERG_FLUXES Compute iceberg fluxes.
%   [QII0,QTI0,QSI0,M0] = GET_ICEBERG_FLUXES(H0,T0,S0,I0,ZI,P)
%   computes the iceberg fluxes for the given parameters.

H0 = s.H(:,i);
T0 = s.T(:,i);
S0 = s.S(:,i);
I0 = s.I(:,i);

if p.M0==0 % if no icebergs

    QIi0 = 0*H0;
    QTi0 = 0*H0;
    QSi0 = 0*H0;
    QVmi0 = 0*H0;
    QTmi0 = 0*H0;
    QSmi0 = 0*H0;
    M0 = 0*H0;

else

    % melting in boxes
    zj = cumsum(H0)-H0/2; % mean depth of boxes
    Tf = p.l1*S0 + p.l2 + p.l3*zj; % local freezing point
    % melt flux
    meltflux = max(0,p.M0*(T0-Tf).*I0);
    meltflux(I0==0) = 0;

    % mixing fluxes between boxes to account for buoyant rising of
    % freshwater
    % Preallocate variables
%     [QVmI, QTmI, QSmI, gmelt, Heff, SA_ice] =  deal(zeros(1, length(I0)-1));
    % effective temperature of meltwater
    Tmelt = -p.l/p.cw;
    % buoyancy difference between melt and ambient
    gmelt = p.g*(p.betaS*S0-p.betaT*(T0-Tmelt));
    % potential upwelling flux
    QVmI = p.U0*p.alphaI^(2/3)*meltflux.^(1/3).*gmelt.^(1/3).*H0.*I0.^(2/3);
   % QVmI(meltflux==0) = 0;
    % scale for density stratification
    gk = max(0,[NaN;p.g*(p.betaS*(S0(2:end)-S0(1:end-1))-p.betaT*(T0(2:end)-T0(1:end-1)))]);
    lengthfac = (1/p.alphaI^(2/3))*((meltflux./I0).^2./(gmelt.*H0)).^(1/3).*gmelt./gk;
    scalefac = 1-exp(-lengthfac);
    scalefac(gk==0) = 1;
    scalefac(I0==0) = 0;
    scalefac(1) = 0; % no upwelling to atmosphere
    QVmI = scalefac.*QVmI;
    % associated temperature and salt fluxes
    QTmI = QVmI.*T0;
    QSmI = QVmI.*S0;

%     % loop over interfaces, allow plumes to rise from below sill layer
%     for k=1:p.N+p.sill-1
%         % reduced gravity
%         gmelt(k) = p.g*(p.betaS*S0(k)-p.betaT*T0(k)); % relative to the iceberg
%         % effective height rise of plume in layer
%         inds = find(zi0<=ints(k) & zi0>=ints(k+1));
%         Heff(k) = trapz(zi0(inds),Ii0(inds));
%         % Surface area of icebergs to melt in layer k (assume characteristic length scale based on Heff(k))
%         % if dynamic, allow MITgcm profile, keep width divided by Heff=
%         % height of box
%         % if static, change to cubes
%         % plus discrete iceberg concentration 
%         SA_ice(k) = 6*sqrt(6)*p.W*p.L*Heff(k)/(sqrt(5*Heff(k)^2/4));
%         % fluxes
%         if Heff(k) == 0 % no iceberg concentration (e.g. on first timestep)
%             QVmI(k) = 0;
%             QTmI(k) = 0;
%             QSmI(k) = 0;
%         elseif gmelt(k) < 0 % iceberg plume is denser than the box
%             QVmI(k) = 0;
%             QTmI(k) = 0;
%             QSmI(k) = 0;
%         else % iceberg plume which mixes into box above 
%             QVmI(k) = (QIi0(k+1))^(1/3)*p.alphaI^(2/3)*gmelt(k)^(1/3)*(SA_ice(k)/Heff(k))^(2/3)*Heff(k)/2;
%             QTmI(k) = QVmI(k)*T0(k+1)+QTi0(k+1);
%             QSmI(k) = QVmI(k)*S0(k+1)+QSi0(k+1);
%         end
% 
%     end

    % final upwelling fluxes
    QVmi0 = [QVmI(2:end)',0]'-[0,QVmI(2:end)']';
    QTmi0 = [QTmI(2:end)',0]'-[0,QTmI(2:end)']';
    QSmi0 = [QSmI(2:end)',0]'-[0,QSmI(2:end)']';

    % assume that the same fraction of meltwater upwells
    QIi0 = (1-scalefac).*meltflux + [scalefac(2:end);0].*[meltflux(2:end);0];
    % associated heat and salt fluxes
    QTi0 = -QIi0*p.l/p.cw;
    QSi0 = -QIi0.*S0;

%     fraction of meltwater stays in box, while fraction rises up as plume
%     QIi0 = (1-p.gamma)*QIi0;
%     QTi0 = (1-p.gamma)*QTi0;
%     QSi0 = (1-p.gamma)*QSi0;
%     QVmi0 = p.gamma*QVmi0;
%     QTmi0 = p.gamma*QTmi0;
%     QSmi0 = p.gamma*QSmi0;

    % 
    M0 = 0*H0;



    % if p.sill % ensures there is no mixing at the sill layer interface
    %     QVmi0(p.N+p.sill) = 0;
    %     QTmi0(p.N+p.sill) = 0;
    %     QSmi0(p.N+p.sill) = 0;
    % end


end

Qi.I = QIi0;
Qi.T = QTi0;
Qi.S = QSi0;
Qi.M = M0;
Qi.Vm = QVmi0;
Qi.Tm = QTmi0;
Qi.Sm = QSmi0;

end