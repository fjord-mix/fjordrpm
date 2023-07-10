function H = get_fjord_boxes_from_density(Tz,Sz,z,p)
    
    inds_interfaces=NaN([p.N-1, 1]);

    %% Define box depths from sigma(T,S) time-averaged profile
    salt_profile=mean(Sz(:,1),2,'omitnan');
    temp_profile=mean(Tz(:,1),2,'omitnan');    

    % will ignore the very top, ensuring minimum layer thickness for box1
    search_T_profile = temp_profile(z < -p.Hmin);
    search_S_profile = salt_profile(z < -p.Hmin);
    search_Z_profile = z(z < -p.Hmin);

    % depth in m is roughly equivalent to pressure in dbar    
    sigma_profile = gsw_sigma0(search_S_profile,search_T_profile);
    
    if ~isfield(p,'sigma_bnds')        
        disp('Warning: sigma values for box interfaces nor defined')
        disp('Using salinity minimum for first interface, sill depth for last')
        disp('any interfaces in-between will be equally spaced')

        if p.N > 1 % there are at least two layers (apart from the possible one below the sill)
            [~,inds_interfaces(1)] = min(search_S_profile); % define one interface as the PW base, regarded as the S minimum        

            % divides the remainder of the water column in equal parts
            z_interface_sfc = search_Z_profile(inds_interfaces(1));
            z_to_divide=p.silldepth-z_interface_sfc;
            delta_z = z_to_divide/(p.N-1); % we subtract 1 because we already have the surface
            for i=2:p.N
                [~,inds_interfaces(i)]=min(abs(search_Z_profile-(z_interface_sfc+(i-1)*delta_z))); % finds closest value
            end
        end
    else
        % define interfaces from sigma intervals
        if p.N-1 > length(p.sigma_bnds)
            disp('Layer interfaces being defined from sigma values!')
            disp('One value per interface (excl. the sill) is needed!')
        end
        for i=1:p.N-1
            [~,inds_interfaces(i)] = min(abs(sigma_profile-p.sigma_bnds(i))); % finds closest value
        end        
    end
    
    z_interfaces=search_Z_profile(inds_interfaces);

    % adds the sill depth if it exists
    if p.sill
        z_interfaces(end+1) = p.silldepth; 
    end
    
    H=z_interfaces(2:end)-z_interfaces(1:end-1);
    H = abs([z_interfaces(1), H,-p.H-z_interfaces(end)]); % adds the bottom of the fjord
    
    for k=2:length(H)-1               % goes down the boxes        
        if H(k) < p.Hmin              % if k-th box is too thin
            h_needed = p.Hmin-H(k);   % checks how much it is missing
            H(k) = H(k)+h_needed;     % adds it up
            H(k+1) = H(k+1)-h_needed; % removes from the box below
            % last box should be substantially thick that no problem would come up
            % However, if that happens, we set it to minimum thickness and
            % display a warning
            if H(k+1) < 0
                H(k+1) = p.Hmin;
                disp(['WARNING: bottom box was too thin and was artificially set to ',num2str(p.Hmin),'m!'])
                disp(["This only happens if the sill is only ",num2str(p.Hmin),"m above the fjord's max depth."])
                disp('This means that the fjord is was made slightly deeper than in reality')
            end
        end                            
    end    
    H = double(H);

% Plot to check the box distribution
% ints=[0,-cumsum(H)];
% figure; hold on;
% plot(sigma_profile,search_Z_profile);
% ylabel('depth (m)'); xlabel('sigma (kg/m3)')
% for j=1:size(ints,2)
%     yline(ints(j));
% end
% ylim([ints(end) ints(1)])

end