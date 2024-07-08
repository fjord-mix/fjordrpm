function [QVi0, QTi0, QSi0, QMi0] = get_iceberg_fluxes(i, p, s)

% GET_ICEBERG_FLUXES Compute iceberg fluxes.
%   [QVi0, QTi0, QSi0, QMi0] = GET_ICEBERG_FLUXES(i, p, s)
%   computes the iceberg fluxes QVi0, QTi0, QSi0, QMi0
%   for the given parameters p at timestep i.

% Get tracer variables at timestep i
H0 = s.H; T0 = s.T(:,i); S0 = s.S(:,i); I0 = s.I;

if p.M0==0 
    % If there are no icebergs, the fluxes are zero by default
    [QVi0, QTi0, QSi0, QMi0] = deal(0*H0);
else
    % Compute the melt flux into each box
    zj = cumsum(H0)-H0/2; % Mean depth of boxes
    Tf = p.l1*S0 + p.l2 + p.l3*zj; % Local freezing point
    QMi0 = max(0,p.M0*(T0-Tf).*I0);

    % Compute velocity scale and entrainment into upwelling
    Teff = -p.l/p.cw; % effective temperature of meltwater   
    gmelt = p.g*(p.betaS*S0-p.betaT*(T0-Teff)); % buoyancy difference
    vel = (QMi0.*gmelt.*H0./(p.alphaI*I0)).^(1/3);
    vel(I0==0 | gmelt<=0) = 0;
    Qent = p.alphaI*vel.*I0;    

    % Compute length scale and fraction for upwelling
    gk = max(0,[p.g*(p.betaS*(S0(2:end)-S0(1:end-1)) ...
              -p.betaT*(T0(2:end)-T0(1:end-1)))]);
    lice = (vel(2:end).^2./H0(2:end)).*(H0(1:end-1)+H0(2:end))./gk;
    lice(vel(2:end)==0) = 0;
    fice = [0;min(1,lice./H0(2:end))];

    % Scale the upwelling flux
    Qentscaled = p.U0*fice.*Qent;

    % Net volume, heat and salt fluxes from upwelling
    QVi0 = [Qentscaled(2:end);0]-[0;Qentscaled(2:end)];
    QTi0 = [Qentscaled(2:end).*T0(2:end);0]-[0;Qentscaled(2:end).*T0(2:end)];
    QSi0 = [Qentscaled(2:end).*S0(2:end);0]-[0;Qentscaled(2:end).*S0(2:end)];

    % Add the meltwater contributions to heat and salt fluxes
    meltcont = [fice(2:end);0].*[QMi0(2:end);0] + (1-fice).*QMi0;
    QTi0 = QTi0-meltcont*p.l/p.cw;
    QSi0 = QSi0-meltcont.*S0;

end

end