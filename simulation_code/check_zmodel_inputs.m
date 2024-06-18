function [status, a] = check_zmodel_inputs(p, a)

% CHECK_ZMODEL_INPUTS Check for errors in zmodel inputs.
%   STATUS = CHECK_ZMODEL_INPUTS(P, A) checks the zmodel inputs P and A,
%   AND returns an error status S with default value 0, or 1 if there was
%   an error. Also returns an updated A if any changes to e.g. dimensions
%   were made.

status = 0;

% Check shelf oscillation parameters have been set up correctly.
if (isfield(p,'zd') && p.zd > 0) && (isfield(p,'tw') && p.tw  <= 0)
    disp('Error: must have positive oscillation period if oscillation strength is set.')
    status = 1;
    return
end

% Check inputs are consistent with specified number of layers.
if any([length(a.H0) ~= p.N,length(a.T0) ~= p.N,length(a.S0) ~= p.N])
    disp('Error: initial conditions not consistent with number of layers');
    status = 1; 
    return
end

% Check sum of layer thicknesses is equal to fjord depth.
if abs(sum(a.H0)-p.H) > 1e-10
    disp('Error: initial box thicknesses must sum to fjord depth');
    status = 1; 
    return
end

% Check the dimensions of the initial conditions and transpose if required.
if size(a.H0, 1) < size(a.H0, 2)
    a.H0 = a.H0';
end

end