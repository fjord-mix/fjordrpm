function [T0, S0] = bin_ocean_profiles(Tz,Sz,z,H0,p)

% calculate mean shelf T/S over box model layers;
ints=[0 -cumsum(H0)];
nan_entries = isnan(Tz) | isnan(Sz);
Tz = Tz(~nan_entries);
Sz = Sz(~nan_entries);
z  =  z(~nan_entries);

Ss = flip(Sz);
Ts = flip(Tz);
zs = flip(z);
zs0 = unique(sort([zs,-cumsum(H0)]));

Ss0 = interp1(zs,Ss,zs0,'pchip','extrap');
Ts0 = interp1(zs,Ts,zs0,'pchip','extrap');

[Te0,Se0] = deal(zeros(size(H0)));

for k=1:length(ints)-1
    inds = find(zs0<=ints(k) & zs0>=ints(k+1));
    if length(inds) == 1 % if there is only one data point, no need to average it
        Se0(k) = Ss0(inds);
        Te0(k) = Ts0(inds);
    elseif H0(k) < 1.1*p.Hmin 
        % sometimes a very thin layer with sharp gradients will
        % not yield satisfactory results. even with a high-res profile
        % so we use a simple average instead of numerical integration
        Se0(k) = mean(Ss0(inds));
        Te0(k) = mean(Ts0(inds));
    else
        Se0(k) = trapz(zs0(inds),Ss0(inds))/H0(k);
        Te0(k) = trapz(zs0(inds),Ts0(inds))/H0(k);
    end
end
S0 = Se0';
T0 = Te0';        

% Plot to check interpolation
% figure;
% subplot(1,2,1); hold on
% plot(Ts0,zs0);
% scatter(T0,(ints(1:end-1)+ints(2:end))/2);
% for j=1:size(ints,2)
%     yline(ints(j));
% end
% ylim([ints(end) ints(1)])
% subplot(1,2,2); hold on
% plot(Ss0,zs0);
% scatter(S0,(ints(1:end-1)+ints(2:end))/2);
% for j=1:size(ints,2)
%     yline(ints(j));
% end
% ylim([ints(end) ints(1)])


end