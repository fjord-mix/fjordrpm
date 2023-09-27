function [T0, S0] = bin_ocean_profiles(Tz,Sz,z,H0,~)
% BIN_OCEAN_PROFILES Bins temperature (Tz) and salinity (Sz) profiles to boxmodel layers (H0)
%   BIN_OCEAN_PROFILES(Tz,Sz,Z,H0,p) returns a "box profile" for
%   temperature (T0) and salinity (S0) for the boxes given by H0

% calculate mean shelf T/S over box model layers;
ints=[0 -cumsum(H0)];
nan_entries = isnan(Tz) | isnan(Sz);
Tz = Tz(~nan_entries);
Sz = Sz(~nan_entries);
z  =  z(~nan_entries);

Ss = flip(Sz);
Ts = flip(Tz);
zs = flip(z);
zs0 = unique(sort([0,zs,-cumsum(H0)]));

Ss0 = interp1(zs,Ss,zs0,'pchip','extrap');
Ts0 = interp1(zs,Ts,zs0,'pchip','extrap');

[Te0,Se0] = deal(zeros(size(H0)));

for k=1:length(ints)-1
    inds = find(zs0<=ints(k) & zs0>=ints(k+1));
    Se0(k) = trapz(zs0(inds),Ss0(inds))/H0(k);
    Te0(k) = trapz(zs0(inds),Ts0(inds))/H0(k);
end
S0 = Se0';
T0 = Te0';        

% Plot to check interpolation
% figure(98);
% subplot(1,2,1); hold on; box on
% plot(Ts0,zs0);
% scatter(T0,(ints(1:end-1)+ints(2:end))/2,'filled');
% for j=1:size(ints,2)
%     yline(ints(j));
% end
% ylim([ints(end) ints(1)])
% ylabel('Depth (m)'); xlabel('Temperature (^oC)')
% subplot(1,2,2); hold on; box on
% plot(Ss0,zs0);
% scatter(S0,(ints(1:end-1)+ints(2:end))/2,'filled');
% for j=1:size(ints,2)
%     yline(ints(j));
% end
% ylim([ints(end) ints(1)])
% xlabel('Salinity (-)')

end