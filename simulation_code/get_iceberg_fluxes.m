function [QVi0, QTi0, QSi0, QMi0] = get_iceberg_fluxes(i, p, s)

% GET_ICEBERG_FLUXES Compute iceberg fluxes.
%   [QVi0, QTi0, QSi0, QMi0] = GET_ICEBERG_FLUXES(i, p, s)
%   computes the iceberg fluxes QVi0, QTi0, QSi0, QMi0
%   for the given parameters p at timestep i.

% Get tracer variables at timestep i
H0 = s.H; T0 = s.T(:,i); S0 = s.S(:,i); I0 = s.I;

if p.M0==0 
    % If there are no icebergs, the fluxes are zero by default
    [QVi0, QTi0, QSi0] = deal(0*H0);
else
    % Compute the melt flux into each box
    zj = cumsum(H0)-H0/2; % Mean depth of boxes
    Tf = p.l1*S0 + p.l2 + p.l3*zj; % Local freezing point
    QMi0 = max(0,p.M0*(T0-Tf).*I0);
    % Set the melt flux to zero where I0 is zero to avoid NaN values
    QMi0(I0==0) = 0;

    % Compute mixing flux QVmI between boxes to account for the buoyant
    % rising of freshwater
    Teff = -p.l/p.cw; % effective temperature of meltwater   
    gmelt = p.g*(p.betaS*S0-p.betaT*(T0-Teff)); % buoyancy difference between melt and ambient
    Qent = p.U0*p.alphaI^(2/3)*QMi0.^(1/3).*gmelt.^(1/3).*H0.*I0.^(2/3);

    % Compute scale factor to account for density stratification between
    % boxes
    scalefac = get_upwelling_scalefactor(p, H0, T0, S0, I0, QMi0, gmelt);
  
    % Scale the mixing flux
    Qent = scalefac.*Qent;
    % Compute the associated temperature and salt fluxes
    QTent = Qent.*T0;
    QSent = Qent.*S0;

    % Calculate the final upwelling fluxes
    QVi0 = [Qent(2:end)',0]'-[0,Qent(2:end)']';
    QTi0 = [QTent(2:end)',0]'-[0,QTent(2:end)']';
    QSi0 = [QSent(2:end)',0]'-[0,QSent(2:end)']';

    % Add the meltwater contributions to heat and salt fluxes
    meltcont = [scalefac(2:end);0].*[QMi0(2:end);0] + (1-scalefac).*QMi0;
    QTi0 = QTi0-meltcont*p.l/p.cw;
    QSi0 = QSi0-meltcont.*S0;

end

end