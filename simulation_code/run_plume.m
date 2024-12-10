function [Qent, Qmelt, knb] = run_plume(j, p, kgl, H0, S0, T0, Qsg0);

% orient properties deepest first
Ta = flipud(T0);
Sa = flipud(S0);
H0 = flipud(H0);
kgl = length(H0)-kgl+1;

% initialise output variables
Qent = zeros(length(H0),1);
Qmelt = zeros(length(H0),1);

% initialise plume model variables
Tp(kgl) = p.l2+p.l3*p.Hgl(j);
Sp(kgl) = 0;
gp(kgl) = p.g*(p.betaS*(Sa(kgl)-Sp(kgl))-p.betaT*(Ta(kgl)-Tp(kgl)));
b(kgl) = (p.alphap*(Qsg0/p.Wp(j))^2/gp(kgl))^(1/3);
u(kgl) = Qsg0/(p.Wp(j)*b(kgl));
edot(kgl) = p.alphap*u(kgl);
QV(kgl) = u(kgl)*b(kgl);
QM(kgl) = u(kgl)^2*b(kgl);
QT(kgl) = b(kgl)*u(kgl)*Tp(kgl);
QS(kgl) = b(kgl)*u(kgl)*Sp(kgl);
[mdot(kgl),Tb(kgl),Sb(kgl)] = meltrate(p,u(kgl),Tp(kgl),Sp(kgl),p.Hgl(j));

% loop over layers
k = kgl;
while gp(k)>0 & k<length(H0)

    % advance the fluxes
    k = k+1;
    if k==kgl+1
        dz = p.Hgl(j)-(p.H-sum(H0(1:kgl)));
    else
        dz = H0(k-1);
    end
    Qent(k-1) = dz*edot(k-1);
    Qmelt(k-1) = dz*mdot(k-1);
    QV(k) = QV(k-1) + Qent(k-1) + Qmelt(k-1);
    QM(k) = QM(k-1) + dz*(b(k-1)*gp(k-1)-p.Cd*u(k-1)^2);
    QT(k) = QT(k-1) + dz*(edot(k-1)*Ta(k-1)+mdot(k-1)*Tb(k-1)-p.Cd^0.5*p.GT*u(k-1)*(Tp(k-1)-Tb(k-1)));
    QS(k) = QS(k-1) + dz*(edot(k-1)*Sa(k-1));

    % update the other required quantities
	b(k) = QV(k)^2/QM(k);
	u(k) = QM(k)/QV(k);
	Tp(k) = QT(k)/QV(k);
	Sp(k) = QS(k)/QV(k);
    gp(k) = p.g*(p.betaS*(Sa(k)-Sp(k))-p.betaT*(Ta(k)-Tp(k)));
    edot(k) = p.alphap*u(k);
    [mdot(k),Tb(k),Sb(k)] = meltrate(p,u(k),Tp(k),Sp(k),p.Hgl(j)-sum(H0(kgl:k-1)));

end

% scale vol fluxes for plume width and orient shallowest first
Qent = flipud(Qent)*p.Wp(j);
Qmelt = flipud(Qmelt)*p.Wp(j);

% properties at neutral buoyancy
knb = length(H0)-length(gp)+1;

end