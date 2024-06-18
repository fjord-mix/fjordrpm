function [QIi0, QTi0, QSi0, QVmi0, QTmi0, QSmi0, M0] = get_zmodel_iceberg_fluxes(i, p, s)

% GET_ZMODEL_ICEBERG_FLUXES Compute mixing fluxes for the zmodel.
%   [QII0, QTI0, QSI0, QVMI0, QTMI0, QSMI0, M0] = GET_ZMODEL_ICEBERG_FLUXES(I, P, S) computes the
%   iceberg fluxes QII0, QTI0, QSI0, QVMI0, QTMI0, QSMI0 and iceberg parameter M0 for the given parameters P and boundary
%   conditions F at timestep I.

% Get tracer variables at timestep i.
H0 = s.H(:,i); T0 = s.T(:,i); S0 = s.S(:,i); I0 = s.I(:,i);

if p.M0==0 
    % If there are no icebergs, the fluxes are zero by default.
    [QIi0, QTi0, QSi0, QVmi0, QTmi0, QSmi0, M0] = deal(0*H0);
else
    % Initialise variables

    
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
%     QVmI(meltflux==0) = 0;
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

    % final upwelling fluxes
    QVmi0 = [QVmI(2:end)',0]'-[0,QVmI(2:end)']';
    QTmi0 = [QTmI(2:end)',0]'-[0,QTmI(2:end)']';
    QSmi0 = [QSmI(2:end)',0]'-[0,QSmI(2:end)']';

    % assume that the same fraction of meltwater upwells
    QIi0 = (1-scalefac).*meltflux + [scalefac(2:end);0].*[meltflux(2:end);0];
    % associated heat and salt fluxes
    QTi0 = -QIi0*p.l/p.cw;
    QSi0 = -QIi0.*S0;

    M0 = 0*H0;

end

end