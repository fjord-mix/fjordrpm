function status = check_zmodel_initialisation(p, a)

status = 0;

% Check shelf oscillation parameters have been set up correctly.
if (isfield(p,'zd') && p.zd > 0) && (isfield(p,'tw') && p.tw  <= 0)
    disp('Error: must have positive oscillation period if oscillation strength is set.')
    status = 1; % Record an error and return error status.
    return
end

% Check initialisation is consistent with specified number of layers.
if any([length(a.H0) ~= p.N,length(a.T0) ~= p.N,length(a.S0) ~= p.N])
    disp('Error: Initial conditions not consistent with number of layers');
    status = 1; 
    return
end

% Check sum of layer thicknesses is equal to fjord depth.
if abs(sum(a.H0)-p.H) > 1e-10
    disp('Error: box thicknesses must sum to fjord depth');
    status = 1; 
    return
end

end
