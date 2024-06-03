function s = check_zmodel(p, H)

    % Check sum of layer thicknesses is equal to fjord depth.
    if abs(sum(H)-p.H) > 1e-10
        disp('Error: box thicknesses must sum to fjord depth');
        s.status = 1; % status == 1 means there was an error
    end

end