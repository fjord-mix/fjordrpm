function a = get_initial_conditions(p, f)

% GET_INITIAL_CONDITIONS Obtain the initial conditions.
%   A = GET_FORCING(P, F) gets the initial conditions A for the fjord for
%   the given box model parameters P and forcing F.

%% Initial fjord layer thicknesses
if p.sill
    if p.fixedthickness==1
        a.H0 = (p.H/(p.N))*ones(1,p.N+1); % no extra layers with a sill 
    else
        a.H0 = [(abs(p.silldepth)/p.N)*ones(1,p.N),p.H-abs(p.silldepth)];
    end
else
    a.H0 = (p.H/p.N)*ones(1,p.N);
end

%% No icebergs in the fjord initially if dynamic, 
% a.I0 = p.icestatic*p.B*f.xi; % 
a.I0 = 0*f.xi;

%% Set initial fjord T/S to be in equilibrium with shelf
ints = [0;-cumsum(a.H0')];
zs0 = unique(sort([f.zs,-cumsum(a.H0)]));
Ss0 = interp1(f.zs,f.Ss(:,1),zs0,'pchip','extrap');
Ts0 = interp1(f.zs,f.Ts(:,1),zs0,'pchip','extrap');
for k=1:length(ints)-1
    inds = find(zs0<=ints(k) & zs0>=ints(k+1));
    a.S0(k) = trapz(zs0(inds),Ss0(inds))/a.H0(k);
    a.T0(k) = trapz(zs0(inds),Ts0(inds))/a.H0(k);
end

end
