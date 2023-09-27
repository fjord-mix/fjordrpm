function [QVb0,QTb0,QSb0] = get_artificial_fluxes(QV0,H0,V0,T0,S0,zs,~,Ss,p)

% GET_ARTIFICIAL_FLUXES Compute artificial fluxes.
%   [QVB0,QTB0,QSB0] = GET_ARTIFICIAL_FLUXES(QV0,H0,V0,T0,S0,ZS,TS,SS,P)
%   computes the artifcial fluxes for the given parameters. 

if p.sill == 1
    QVsill = QV0(p.N+p.sill);
    QTsill = (QVsill>0)*QVsill*T0(p.N+p.sill)+(QVsill<0)*QVsill*T0(p.N);
    QSsill = (QVsill>0)*QVsill*S0(p.N+p.sill)+(QVsill<0)*QVsill*S0(p.N);
    % resulting box-to-box vectors
    QVb0 = [zeros(p.N-p.sill,1);QVsill;-QVsill];
    QTb0 = [zeros(p.N-p.sill,1);QTsill;-QTsill];
    QSb0 = [zeros(p.N-p.sill,1);QSsill;-QSsill];
else
    QVb0 = 0*H0;
    QTb0 = 0*H0;
    QSb0 = 0*H0;
end

%% do layer nudging followed by minimum thickness
if ~isnan(p.trelax) % if layer nudging active
    [QVb0,QTb0,QSb0]=get_nudging_fluxes(QVb0,QTb0,QSb0,H0,T0,S0,Ss,zs,p);
end

if ~isnan(p.Hmin) % if minimum thickness active
    [QVb0, QTb0, QSb0] = get_hmin_fluxes(QV0,QVb0,QTb0,QSb0,V0,H0,T0,S0,p);
end

end