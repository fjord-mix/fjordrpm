function [QVs0,QTs0,QSs0,Se0,Te0,phi0] = get_zmodel_shelf_fluxes(i, p, f, s)

% GET_SHELF_FLUXES Compute shelf fluxes.
%   [QVS0,QTS0,QSS0,SE0,TE0,PHI0] = GET_SHELF_FLUXES(H0,T0,S0,ZS,TS,SS,QSG0,P)
%   computes the shelf fluxes for the given parameters. 

H0 = s.H(:,i);
T0 = s.T(:,i);
S0 = s.S(:,i);
Qsg0 = f.Qsg(i);
zs = f.zs;
Ts = f.Ts(:,i);
Ss = f.Ss(:,i);


if p.C0==0
    QVs0 = 0*H0;
    QTs0 = 0*H0;
    QSs0 = 0*H0;
    Se0 = NaN*H0;
    Te0 = NaN*H0;
    phi0 = NaN*H0;

else
    % calculate mean shelf T/S over box model layers above sill
    [gp, phi0] = deal(zeros(p.N, 1));
    [Te0, Se0] = bin_ocean_profiles(Ts,Ss,zs,H0',p);

    % get fjord to shelf reduced gravity
    for k=1:s.ksill
        gp(k) = p.g*(p.betaS*(S0(k)-Se0(k))-p.betaT*(T0(k)-Te0(k)));
    end

    % calculate potentials over above-sill layers
    for k=1:s.ksill
        if k==1
            phi0(k) = gp(k)*H0(k)/2;
        else
            phi0(k) = phi0(k-1)+gp(k-1)*H0(k-1)/2+gp(k)*H0(k)/2;
        end
    end

    % above-sill fluxes before barotropic compensation
    Q = p.C0*p.W*H0(1:s.ksill).*phi0(1:s.ksill)/p.L;

%     if p.fixedthickness==1 && p.sill==1
%         % no exchange below the sill height
%         Q(-cumsum(H0)< -p.silldepth) = 0;
%     end

    % above-sill fluxes after ensuring depth mean = QSg0 when plume is turned on,
    % and 0 when plume is turned off
    if p.P0==0
        Qsg0 = 0;
    end
    % above-sill fluxes after barotropic compensation
    QVs0 = Q + H0(1:s.ksill)*(Qsg0-sum(Q))/sum(H0(1:s.ksill));
%     if p.sill==1
%         if p.fixedthickness==0
%             QVs0(p.N+p.sill) = 0;
%             Se0(p.N+p.sill) = 0;
%             Te0(p.N+p.sill) = 0;
%         elseif p.fixedthickness==1
%             % only computed based on the above-sill layers
%             QVs0 = Q + [H0(-cumsum(H0)>=-p.silldepth); 0*H0(-cumsum(H0)<-p.silldepth)]*(Qsg0-sum(Q))/sum([H0(-cumsum(H0)>=-p.silldepth); 0*H0(-cumsum(H0)<-p.silldepth)]);
%         end
%     end
    % fill any below-sill exchanges with zeros
    QVs0(s.ksill+1:p.N) = 0;
   
    % resulting heat/salt fluxes
    QTs0 = (QVs0>=0).*QVs0.*T0 + (QVs0<0).*QVs0.*Te0;
    QSs0 = (QVs0>=0).*QVs0.*S0 + (QVs0<0).*QVs0.*Se0;
end

Qs.V = QVs0;
Qs.T = QTs0;
Qs.S = QSs0;
Qs.Se = Se0;
Qs.Te = Te0;
Qs.phi = phi0;

end