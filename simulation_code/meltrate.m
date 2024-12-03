function [mdot,Tb,Sb] = meltrate(p, u0, T0, S0, z0)

% meltrate Compute submarine melt rate.
%   mdot = meltrate(p, T0, S0, z0) computes the submarine melt rate (m/s)
%   using the 3-equation parameterisation.

a1 = p.l1*(p.ci*p.GS-p.cw*p.GT);
a2 = p.cw*p.GT*(T0-p.l2-p.l3*z0)+p.GS*(p.ci*(p.l2+p.l3*z0-p.l1*S0-p.Ti)+p.l);
a3 = -p.GS*S0.*(p.ci*(p.l2+p.l3*z0-p.Ti)+p.l);
Sb = (-a2+sqrt(a2.^2-4*a1*a3))./(2*a1);
Tb = p.l1*Sb+p.l2+p.l3*z0;
mdot = p.cw*p.Cd^(1/2)*p.GT*u0.*(T0-Tb)./(p.l+p.ci*(Tb-p.Ti));

end