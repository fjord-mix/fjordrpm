function a = get_initial_conditions(p, f)

% GET_INITIAL_CONDITIONS Obtain the initial conditions.
%   A = GET_INITIAL_CONDITIONS(P, F) gets the initial conditions A for the fjord for
%   the given box model parameters P and forcing F.

%% Set the initial fjord layer thicknesses.
a.H0 = (p.H/p.N)*ones(1,p.N);


%% Initial icebergs in fjord
if p.icestatic
    % a.I0 is the surface area of icebergs in a box
    % so either use the idealised expression here or load from file
    a.I0 = p.A0*p.if(p.nu0, abs(p.zgl), -cumsum(a.H0)+a.H0/2);    
else
    a.I0 = 0*a.H0;
end



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
