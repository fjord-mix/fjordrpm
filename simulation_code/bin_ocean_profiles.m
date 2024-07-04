function [Te0, Se0] = bin_ocean_profiles(Tz, Sz, z, H0)

% BIN_OCEAN_PROFILES Interpolates shelf profile to layers.
%   [Te0, Se0] = BIN_OCEAN_PROFILES(Tz, Sz, z, H0) bins temperature (Tz)
%   and salinity (Sz) profiles to model layers (H0) and returns a "layer
%   profile" for temperature (Te0) and salinity (Se0).

% Remove any nan entries from the shelf temperature, salinity and depth
% profiles
nan_entries = isnan(Tz) | isnan(Sz);
Ts = Tz(~nan_entries);
Ss = Sz(~nan_entries);
zs = z(~nan_entries);

% Add the layer boundaries H0 into the vector of shelf z-values zs
zs0 = unique(sort([0; zs'; -cumsum(H0)]));

% Interpolate the shelf temperature and salinitity profiles given on sample
% points zs onto the new grid zs0, allowing for the new grid to expand
% beyond the domain of zs
Ss0 = interp1(zs,Ss,zs0,'pchip','extrap');
Ts0 = interp1(zs,Ts,zs0,'pchip','extrap');

% Calculate mean shelf temperature and salinity Te0 and Se0 over zmodel
% layers
ints=[0; -cumsum(H0)];
[Te0, Se0] = deal(zeros(size(H0)));
for k=1:length(ints)-1
    % Find the boundaries of the box in the shelf grid zs0
    inds = find(zs0<=ints(k) & zs0>=ints(k+1));
    % Average the temperature and salinity profiles over this box
    Se0(k) = trapz(zs0(inds),Ss0(inds))/H0(k);
    Te0(k) = trapz(zs0(inds),Ts0(inds))/H0(k);
end
      
end