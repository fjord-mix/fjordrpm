function [T0, S0, He0] = bin_ocean_profiles_withsill(Tz,Sz,z,H0,p)
% BIN_OCEAN_PROFILES Bins temperature (Tz) and salinity (Sz) profiles to boxmodel layers (H0)
% takes into account the sill, with separate calculations for partially-
% occluded and below-sill layers.
%   BIN_OCEAN_PROFILES(Tz,Sz,Z,H0,p) returns a "box profile" for
%   temperature (T0) and salinity (S0) for the boxes given by H0 and the
%   sill depth

% Calculate mean shelf T/S over box model layers;
ints=[0 -cumsum(H0)]; % First entry surface, last entry bottom depth
nan_entries = isnan(Tz) | isnan(Sz); % get rid of nan entries 
Tz = Tz(~nan_entries);
Sz = Sz(~nan_entries);
z  =  z(~nan_entries);

Ss = flip(Sz); % first entry is surface, last entry is bottom depth
Ts = flip(Tz);
zs = flip(z); % first entry is surface, last entry is bottom depth
zs0 = unique(sort([0,zs,-cumsum(H0)])); % first entry is bottom depth, last entry is surface

% interpolate the shelf data to the box model depth values
Ss0 = interp1(zs,Ss,zs0,'pchip','extrap'); % first entry is bottom depth, last entry is surface 
Ts0 = interp1(zs,Ts,zs0,'pchip','extrap');

% H0 is location of box model layer boundaries 
[Te0,Se0, He0] = deal(zeros(size(H0)));

for k=1:length(ints)-1
    
    % find the corresponding shelf conditions
    % H0 has entries going surface to bottom

    % if the layer is above the sill, integrate as normal
    if abs(ints(k)) < p.silldepth && abs(ints(k+1)) <= p.silldepth
        % above sill
        % integral of the the layer depth on the sill, over the depth of the layer
        % find z bounding values of the layers, first loop is for top layer beneath surface
        inds = find(zs0<=ints(k) & zs0>=ints(k+1));
        Se0(k) = trapz(zs0(inds),Ss0(inds))/H0(k);       
        Te0(k) = trapz(zs0(inds),Ts0(inds))/H0(k);
        He0(k) = H0(k);

    % if the layer is partially occluded, integrate only over the part in
    % contact with the fjord layer (above the sill)
    elseif abs(ints(k)) < p.silldepth && abs(ints(k+1)) > p.silldepth
        % partially occluded
        % z bounding values are only to the top of the sill
        inds = find(zs0 <=ints(k) & zs0>= -p.silldepth);
        % integral of the part of the layer in contact with the fjord, over the depth of the layer
        Se0(k) = trapz(zs0(inds),Ss0(inds))/abs(zs0(inds(1))-zs0(inds(end)));       
        Te0(k) = trapz(zs0(inds),Ts0(inds))/abs(zs0(inds(1))-zs0(inds(end)));
        He0(k) = abs(zs0(inds(1))-zs0(inds(end)));

    % if the layer is completely below the sill, use the reference point on
    % the shelf immediately above the sill
    else 
        % below sill
        % z bounding values are the reference point on the shelf
        % immediately above the sill
        inds = find(zs0== -p.silldepth);
        % No integral- just use reference point salinity 
        Se0(k) = Ss0(inds);       
        Te0(k) = Ts0(inds);
        He0(k) = He0(k-1)/10; % small gap in layer above for flow to go through - check with Donald 
    

    end
% First entry of se0, te0 is for surface layer, last entry is bottom depth 


end
S0 = Se0';
T0 = Te0';
He0 = He0';

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