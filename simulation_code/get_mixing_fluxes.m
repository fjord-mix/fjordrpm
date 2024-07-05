function [QVk0,QTk0,QSk0] = get_mixing_fluxes(i, p, s)

% GET_MIXING_FLUXES Compute vertical tracer mixing fluxes.
%   [QVK0, QTK0, QSK0] = GET_MIXING_FLUXES(i, p, s) computes the
%   vertical mixing fluxes QVK0, QTK0, QSSK0 for the given parameters p
%   and solution s at timestep i.

% Required variables at timestep i
H0 = s.H;
T0 = s.T(:,i);
S0 = s.S(:,i);

% The net volume mixing fluxes are always zero
QVk0 = 0*H0;

% Reduced gravity between layers
gp = p.g*(p.betaS*(S0(2:end)-S0(1:end-1))-p.betaT*(T0(2:end)-T0(1:end-1)));

% Horizontal velocity of layers
u = (s.QVp(:,i)-s.QVs(:,i))./(2*p.W*H0);

% Richardson number
du = u(2:end)-u(1:end-1);
Ri = gp.*(H0(2:end)+H0(1:end-1))./(2*du.^2);
Ri(du==0) = p.Ri0;
Ri(Ri>p.Ri0) = p.Ri0;
Ri(Ri<0) = 0;

% Get diffusivity as a function of the Richardson number
Kz = (Ri<p.Ri0 & Ri>0)*p.K0.*(1-(Ri/p.Ri0).^2).^3;

% Compute the mixing fluxes going in/out of each layer.
QS = 2*p.W*p.L*Kz.*(S0(2:end)-S0(1:end-1))./(H0(2:end)+H0(1:end-1));
QT = 2*p.W*p.L*Kz.*(T0(2:end)-T0(1:end-1))./(H0(2:end)+H0(1:end-1));

% The final layer fluxes are the net of the interface fluxes
QTk0 = [QT;0]-[0;QT];
QSk0 = [QS;0]-[0;QS];

end