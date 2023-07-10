function Sint =  get_interface_salinities(zi,Ti,Si,p)
    % getting current "shelf box profile" for nudging  
    % z_interfaces=-cumsum(Hi(1:end-1-has_sill));    
    % Sint = interp1(zi,Si,z_interfaces,"linear");    

    % sill_depth=-cumsum(Hi(1:end-1-has_sill));
    H_shelf = double(get_fjord_boxes_from_density(Ti,Si,zi,p));
    [~,S_shelf] = bin_ocean_profiles(Ti,Si,zi,H_shelf,p);
    % Sint = S_shelf;%(1:end-1);
    Sint = 0.5*(S_shelf(2:end)+S_shelf(1:end-1));
end