function status = check_inputs(p, t, f, a)

% CHECK_INPUTS Check for errors in model inputs.
%   status = CHECK_INPUTS(p, t, f, a) checks the model inputs.

status = 0;

% Check provided depths are positive
if any([p.Hsill<0, p.Hgl<0, p.H<0])
    disp('Error: p.H, p.Hgl and p.Hsill must be positive');
    status = 1;
    return
end

% Check dimensionality of initial conditions
if ~isequal([size(a.H0)],[size(a.S0)],[size(a.T0)],[p.N,1])
    disp('Error: initial conditions must have dimensions p.N x 1');
    status = 1;
    return
end

% Check sum of layer thicknesses is equal to fjord depth
if abs(sum(a.H0)-p.H) > 1e-10
    disp('Error: layer thicknesses (a.H0) must sum to fjord depth (p.H)');
    status = 1; 
    return
end

% Check dimensionality of forcings
nz = length(f.zs);
nt = length(t);
% f.zs
if ~isequal(size(f.zs),[nz,1])
    disp('Error: f.zs must have dimensions nz x 1');
    status = 1;
    return
end
% f.Ss and f.Ts
if ~isequal(size(f.Ss),size(f.Ts),[nz,nt])
    disp('Error: f.Ss and f.Ts must have dimensions nz x nt');
    status = 1;
    return
end
% f.Qsg
if ~isequal(size(f.Qsg),[1,nt])
    disp('Error: f.Qsg must have dimensions 1 x nt');
    status = 1;
    return
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

end

