function status = check_zmodel(p, H)

% CHECK_ZMODEL Check for errors in the layer depths.
%   STATUS = CHECK_ZMODEL(I, P, S) checks the zmodel layers H at timestep
%   parameters P and returns an error status S with default value 0, or 1
%   if there was an error.

status = 0;

% Check sum of layer thicknesses is equal to fjord depth.
if abs(sum(H)-p.H) > 1e-10
    disp('Error: box thicknesses must sum to fjord depth');
    status = 1; %
end

% Check that the layer depths are not changing size.
if numel(uniquetol(H,1e-8)) > 1+p.sill
    disp("Layer depths not all the same")
    status = 1;
end

end