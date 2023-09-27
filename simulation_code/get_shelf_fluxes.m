function [QVs0,QTs0,QSs0,Se0,Te0,phi0] = get_shelf_fluxes(H0,T0,S0,zs,Ts,Ss,Qsg0,p)

% GET_SHELF_FLUXES Compute shelf fluxes.
%   [QVS0,QTS0,QSS0,SE0,TE0,PHI0] = GET_SHELF_FLUXES(H0,T0,S0,ZS,TS,SS,QSG0,P)
%   computes the shelf fluxes for the given parameters. 

if p.C0==0
    QVs0 = 0*H0;
    QTs0 = 0*H0;
    QSs0 = 0*H0;
    Se0 = NaN*H0;
    Te0 = NaN*H0;
    phi0 = NaN*H0;

else
    % calculate mean shelf T/S over box model layers above sill
    [gp, phi0] = deal(zeros(1, p.N));
    [Te0, Se0] = bin_ocean_profiles(Ts,Ss,zs,H0',p);

    % get fjord to shelf reduced gravity
    for k=1:p.N
        gp(k) = p.g*(p.betaS*(S0(k)-Se0(k))-p.betaT*(T0(k)-Te0(k)));
    end

    % calculate potentials over above-sill layers
    for k=1:p.N
        if k==1
            phi0(k) = gp(k)*H0(k)/2;
        else
            phi0(k) = phi0(k-1)+gp(k-1)*H0(k-1)/2+gp(k)*H0(k)/2;
        end
    end

    % fluxes before barotropic compensation
    Q = p.C0*p.W*H0(1:p.N).*phi0'/p.L;

    % fluxes after ensuring depth mean = QSg0 when plume is turned on,
    % and 0 when plume is turned off
    if p.P0==0
        Qsg0 = 0;
    end
    QVs0 = Q + H0(1:p.N)*(Qsg0-sum(Q))/sum(H0(1:p.N));
    if p.sill==1
        QVs0(p.N+p.sill) = 0;
        Se0(p.N+p.sill) = 0;
        Te0(p.N+p.sill) = 0;
    end

    % resulting heat/salt fluxes
    QTs0 = (QVs0>0).*QVs0.*T0 + (QVs0<0).*QVs0.*Te0;
    QSs0 = (QVs0>0).*QVs0.*S0 + (QVs0<0).*QVs0.*Se0;
end

end