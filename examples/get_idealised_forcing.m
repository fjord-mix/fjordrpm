function f = get_idealised_forcing(p, t)

% GET_IDEALISED_FORCING  Obtain idealised boundary conditions in functional
% form.
%   F = GET_IDEALISED_FORCING(P, T) gets the forcing structure F for the
%   shelf, icebergs and subglacial discharge for the given box model
%   parameters P and T.

%% Get the functional form of the shelf forcing.
f.zs = linspace(-p.H, 0, abs(p.H)+1); 
if p.zd == 0
    f.zo = p.z0*ones(size(t)); % no oscillation
else
    f.zo = p.z0+(p.zd/2)*sin(2*pi*t/p.tw); % oscillation
end
[ZI, ZS] = meshgrid(f.zo, f.zs);
f.Ss = p.sf(p.Sbottom, p.Stop, ZS, ZI);
f.Ts = p.sf(p.Tbottom, p.Ttop, ZS, ZI); 

%% Get the functional form of the iceberg forcing.
f.zi = f.zs; 
% f.D = p.df(t, p.D0);
f.D = 0*t;
f.xi = p.if(p.nu0, abs(p.zgl), f.zi)';

%% Get the functional form of the subglacial discharge forcing.
f.Qsg = p.gf(t, p.Qv0);

end