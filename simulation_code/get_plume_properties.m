function [kgl, knb, Qp, Sp, Tp] = get_plume_properties(p, H0, S0, T0, Qsg0)

% GET_PLUME_PROPERTIES Compute plume properties for the zmodel.
%   [KGL, KNB, QP, SP, TP] = GET_PLUME_PROPERTIES(P, H0, S0, T0, QSG0)
%   computes the plume fluxes for the given parameters P, zmodel tracers
%   H0, S0, T0 and subglacial discharge flux QSG0.


% Initialise variables
[Qp, Sp, Tp, gp] = deal(zeros(p.N, 1));

% Find the zmodel layer containing the grounding line.
ints = cumsum(H0);
kgl = find(ints >= abs(p.zgl)-1e-6, 1);

% Set the plume properties at the grounding line.
Qp(kgl) = Qsg0;
Sp(kgl) = 0;
Tp(kgl) = 0;
gp(kgl) = p.g*(p.betaS*(S0(kgl)-Sp(kgl))-p.betaT*(T0(kgl)-Tp(kgl)));

% If it exists, the properties at the first interface above grounding
% line need special treatment because the box might be partial if the
% grounding line does not coincide with a box boundary.
if kgl>1
    k = kgl-1;
    Qp(k) = Qp(k+1) + p.P0^(2/3)*Qp(k+1)^(1/3)*gp(k+1)^(1/3)*(abs(p.zgl)-ints(k));
    Tp(k) = (Qp(k+1)*Tp(k+1)+(Qp(k)-Qp(k+1))*T0(k+1))/Qp(k);
    Sp(k) = (Qp(k+1)*Sp(k+1)+(Qp(k)-Qp(k+1))*S0(k+1))/Qp(k);
    gp(k) = p.g*(p.betaS*(S0(k)-Sp(k))-p.betaT*(T0(k)-Tp(k)));
end

% Apply to successive interfaces higher provided the plume is still rising.
while gp(k)>0 && k>1
    k = k-1;
    Qp(k) = Qp(k+1) + p.P0^(2/3)*Qp(k+1)^(1/3)*gp(k+1)^(1/3)*H0(k+1);
    Tp(k) = (Qp(k+1)*Tp(k+1)+(Qp(k)-Qp(k+1))*T0(k+1))/Qp(k);
    Sp(k) = (Qp(k+1)*Sp(k+1)+(Qp(k)-Qp(k+1))*S0(k+1))/Qp(k);
    gp(k) = p.g*(p.betaS*(S0(k)-Sp(k))-p.betaT*(T0(k)-Tp(k)));
end

% Find the neutral buoyancy box, which is set by default to the lowest box.
knb = find(gp<0);
if isempty(knb)
    knb=1;
end

end