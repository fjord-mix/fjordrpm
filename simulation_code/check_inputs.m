function [status, a] = check_inputs(p, t, a)

% CHECK_INPUTS Check for errors in model inputs.
%   [status, a] = CHECK_INPUTS(p, a) checks the zmodel inputs p and a, and
%   returns an error status with default value 0, or 1 if there was an
%   error. Also returns an updated a if any changes were made.

status = 0;

% Check provided depths are positive
if any([p.Hsill<0, p.Hgl<0, p.H<0])
    disp('Error: p.H, p.Hgl and p.Hsill must be positive');
    status = 1;
    return
end

% Check inputs are consistent with specified number of layers
if any([length(a.H0) ~= p.N,length(a.T0) ~= p.N,length(a.S0) ~= p.N])
    disp('Error: initial conditions not consistent with number of layers');
    status = 1; 
    return
end

% Check sum of layer thicknesses is equal to fjord depth
if abs(sum(a.H0)-p.H) > 1e-10
    disp('Error: initial box thicknesses must sum to fjord depth');
    status = 1; 
    return
end

% Check the dimensions of the initial conditions and transpose if required
if size(a.H0, 1) < size(a.H0, 2)
    a.H0 = a.H0';
end

% If t_save is not given, it will take t instead
if ~isfield(p,'t_save')
    p.t_save=t;
end

% Check that the time values where the solution will be saved are a subset
% of the time values where the solution will be computed
if ~all(ismember(p.t_save, t)) 
    disp("Error: the values where the solution is saved, p.t_save, must be a subset of the values where the solution is computed, t.")
    status = 1;
    return
end

