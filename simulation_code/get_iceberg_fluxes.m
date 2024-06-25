function [QIi0, QTi0, QSi0, QVmi0, QTmi0, QSmi0] = get_iceberg_fluxes(i, p, s)

% GET_ICEBERG_FLUXES Compute mixing fluxes for the zmodel.
%   [QII0, QTI0, QSI0, QVMI0, QTMI0, QSMI0, M0] = GET_ICEBERG_FLUXES(I, P,
%   S) computes the iceberg fluxes QII0, QTI0, QSI0, QVMI0, QTMI0, QSMI0
%   for the given parameters P and boundary conditions F at timestep I.

% Get tracer variables at timestep i.
H0 = s.H(:,i); T0 = s.T(:,i); S0 = s.S(:,i); I0 = s.I(:,i);

if p.M0==0 
    % If there are no icebergs, the fluxes are zero by default.
    [QIi0, QTi0, QSi0, QVmi0, QTmi0, QSmi0] = deal(0*H0);
else
    % Compute the melt flux into each box.
    zj = cumsum(H0)-H0/2; % Mean depth of boxes.
    Tf = p.l1*S0 + p.l2 + p.l3*zj; % Local freezing point.
    meltflux = max(0,p.M0*(T0-Tf).*I0);
    % Set the melt flux to zero where I0 is zero to avoid NaN values.
    meltflux(I0==0) = 0;

    % Compute mixing flux QVmI between boxes to account for the buoyant
    % rising of freshwater.
    Tmelt = -p.l/p.cw; % effective temperature of meltwater   
    gmelt = p.g*(p.betaS*S0-p.betaT*(T0-Tmelt)); % buoyancy difference between melt and ambient
    QVmI = p.U0*p.alphaI^(2/3)*meltflux.^(1/3).*gmelt.^(1/3).*H0.*I0.^(2/3);

    % Compute scale factor to account for density stratification between
    % boxes.
    scalefac = get_upwelling_scalefactor(p, H0, T0, S0, I0, meltflux, gmelt);
  
    % Scale the mixing flux.
    QVmI = scalefac.*QVmI;
    % Compute the associated temperature and salt fluxes.
    QTmI = QVmI.*T0;
    QSmI = QVmI.*S0;

    % Calculate the final upwelling fluxes.
    QVmi0 = [QVmI(2:end)',0]'-[0,QVmI(2:end)']';
    QTmi0 = [QTmI(2:end)',0]'-[0,QTmI(2:end)']';
    QSmi0 = [QSmI(2:end)',0]'-[0,QSmI(2:end)']';

    % Assume that the same fraction of meltwater upwells.
    QIi0 = (1-scalefac).*meltflux + [scalefac(2:end);0].*[meltflux(2:end);0];
    % Compute the associated heat and salt fluxes
    QTi0 = -QIi0*p.l/p.cw;
    QSi0 = -QIi0.*S0;

end

end