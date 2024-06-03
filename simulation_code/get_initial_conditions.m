function a = get_initial_conditions(p, f)

% GET_INITIAL_CONDITIONS Obtain the initial conditions.
%   A = GET_FORCING(P, F) gets the initial conditions A for the fjord for
%   the given box model parameters P and forcing F.

%% Set the initial fjord layer thicknesses.
a.H0 = (p.H/p.N)*ones(1,p.N);
if p.sill==1
    if p.fixedthickness==1
        % If the layers are fixed thickness, make sure the sill height corresponds
        % with a layer boundary.
        [minValue, closestIndex] = min(abs(cumsum(a.H0)-p.silldepth));
        a.H0(closestIndex) = a.H0(closestIndex) + minValue;
        a.H0(closestIndex+1) = a.H0(closestIndex+1) - minValue;
    else 
        % If the layers are variable thickness, add an extra layer for the
        % sill.
        a.H0 = [(abs(p.silldepth)/p.N)*ones(1,p.N),p.H-abs(p.silldepth)];
    end    
end

%% No icebergs in the fjord initially if dynamic.
% a.I0 = p.icestatic*p.B*f.xi; % 
a.I0 = 0*f.xi;

%% Set the initial fjord T/S to be in equilibrium with shelf.
% Interpolate the boundary conditions to the model z-discretisation.
zs0 = unique(sort([f.zs,-cumsum(a.H0)]));
Ss0 = interp1(f.zs,f.Ss(:,1),zs0,'pchip','extrap');
Ts0 = interp1(f.zs,f.Ts(:,1),zs0,'pchip','extrap');
% Find the layer boundaries.
ints = [0;-cumsum(a.H0')];
% Bin the boundary conditions to the discrete layers.
for k=1:length(ints)-1
    inds = find(zs0<=ints(k) & zs0>=ints(k+1));
    a.S0(k) = trapz(zs0(inds),Ss0(inds))/a.H0(k);
    a.T0(k) = trapz(zs0(inds),Ts0(inds))/a.H0(k);
end

end
