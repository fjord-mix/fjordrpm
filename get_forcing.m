function f = get_forcing(p, t)

% GET_FORCING  Obtain the forcing structure.
%   F = GET_FORCING(P, T) gets the forcing structure F for the shelf, icebergs and
% subglacial discharge for the given box model parameters P and T.

%% Shelf forcing
f.zs = -p.H:0; % will this break if p.H is not an integer?
if p.zd == 0
    f.zi = p.z0*ones(size(t)); % no oscillation
else
    f.zi = p.z0+(p.zd/2)*sin(2*pi*t/p.tw); % oscillation
end
[ZI, ZS] = meshgrid(f.zi, f.zs);
f.Ss = p.sf(p.Sbottom, p.Stop, ZS, ZI);
f.Ts = 0*ZS + 3; % 3 should be a user-set parameter? 

%% Iceberg forcing
f.zi = (-p.H:10:0)'; % 10 should be a user set parameter?
f.D = zeros(1,length(t));
f.xi = (p.nu0/p.H)*exp(p.nu0*f.zi/p.H)/(1-exp(-p.nu0));

%% Subglacial discharge forcing
f.Qsg = zeros(1,length(t));
f.Qsg(t>5) = p.Qv0; % 5 should be a user-set parameter? 

end