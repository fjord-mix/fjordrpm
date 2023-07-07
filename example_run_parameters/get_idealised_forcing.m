function f = get_idealised_forcing(p, t)

% GET_IDEALISED_FORCING  Obtain the forcing structure for an idealised input of functional form.
%   F = GET_IDEALISED_FORCING(P, T) gets the forcing structure F for the shelf, icebergs and
% subglacial discharge for the given box model parameters P and T.

%% Shelf forcing
f.zs = linspace(-p.H, 0, abs(p.H)+1); 
if p.zd == 0
    f.zi = p.z0*ones(size(t)); % no oscillation
else
    f.zi = p.z0+(p.zd/2)*sin(2*pi*t/p.tw); % oscillation
end
[ZI, ZS] = meshgrid(f.zi, f.zs);
f.Ss = p.sf(p.Sbottom, p.Stop, ZS, ZI);
f.Ts = p.sf(p.Tbottom, p.Ttop, ZS, ZI); 

%% Iceberg forcing
f.zi = f.zs; 
f.D = zeros(1,length(t));
f.xi = p.if(p.nu0, p.H, f.zi);

%% Subglacial discharge forcing
f.Qsg = p.gf(t, p.Qv0);

end