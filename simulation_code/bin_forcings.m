function [Ts, Ss, Qsg] = bin_forcings(f, H, t)

% BIN_FORCINGS Puts forcings on model layers and time steps
%   [Ts, Ss, Qsg] = BIN_FORCINGS(f, H, t) calculates the mean value of
% shelf temperature and salinity over the depths of each model layer. It
% then interpolates both the resulting shelf profiles and subglacial
% discharge onto the model time steps.

% Remove any nan entries from the shelf profiles
% (Assumes the nan entries are the same at every time)
nan_entries = isnan(f.Ts(:,1)) | isnan(f.Ss(:,1));
f.Ts = f.Ts(~nan_entries,:);
f.Ss = f.Ss(~nan_entries,:);
f.zs = f.zs(~nan_entries);

%% First put shelf forcing on model layers

% Add the layer boundaries H0 into the vector of shelf z-values
z0 = unique(sort([0; f.zs; -cumsum(H)]));

% Interpolate the shelf temperature and salinity profiles given on sample
% points zs onto the new grid z0
S0 = interp1(f.zs,f.Ss,z0,'pchip','extrap');
T0 = interp1(f.zs,f.Ts,z0,'pchip','extrap');

% Calculate shelf temperature and salinity Ts0 and Ss0 on model layers
ints = [0; -cumsum(H)];
for k=1:length(ints)-1
    % Find the boundaries of the layer in the shelf grid z0
    inds = find(z0<=ints(k) & z0>=ints(k+1));
    % Average the temperature and salinity profiles over this layer
    Ss(k,:) = trapz(z0(inds),S0(inds,:))/H(k);
    Ts(k,:) = trapz(z0(inds),T0(inds,:))/H(k);
end

%% Second put forcings on model time steps

% Shelf conditions
Ss = interp1(f.ts,Ss',t,'linear')';
Ts = interp1(f.ts,Ts',t,'linear')';

% Subglacial discharge
for j=1:size(f.Qsg,1)
    Qsg(j,:) = interp1(f.tsg,f.Qsg(j,:),t,'linear');
end
      
end