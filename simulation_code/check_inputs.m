function status = check_inputs(p, t, f, a)

% CHECK_INPUTS Check for errors in model inputs.
%   status = CHECK_INPUTS(p, t, f, a) checks the model inputs.

status = 0;

% Check provided depths are positive
if any([p.Hsill<0, p.Hgl<0, p.H<0])
    error('p.H, p.Hgl and p.Hsill must be positive');
end

% Check dimensionality of initial conditions
if ~isequal([size(a.H0)],[size(a.S0)],[size(a.T0)],[size(a.I0)],[p.N,1])
    error('Initial conditions (a.H0, a.S0, a.T0, a.I0) must have dimensions p.N x 1');
end

% Check sum of layer thicknesses is equal to fjord depth
if abs(sum(a.H0)-p.H) > 1e-10
    error('Layer thicknesses (a.H0) must sum to fjord depth (p.H)');
end

% Check dimensionality of shelf forcing
nz = length(f.zs);
nt = length(f.ts);
% f.ts
if ~isequal(size(f.ts),[1,nt])
    error('f.ts must have dimensions 1 x nt');
end
% f.zs
if ~isequal(size(f.zs),[nz,1])
    error('f.zs must have dimensions nz x 1');
end
% f.Ss and f.Ts
if ~isequal(size(f.Ss),size(f.Ts),[nz,nt])
    error('f.Ss and f.Ts must have dimensions length(f.zs) x length(f.ts)');
end

% Check dimensionality of dischaarge forcing
nt = length(f.tsg);
% f.tsg
if ~isequal(size(f.tsg),[1,nt])
    error('f.tsg must have dimensions 1 x nt');
end
% f.Qsg
if ~isequal(size(f.Qsg),[1,nt])
    error('f.Qsg must have dimensions 1 x nt');
end

% % If t_save is not given, it will take t instead
% if ~isfield(p,'t_save')
%     p.t_save=t;
% end
% 
% % Check that the time values where the solution will be saved are a subset
% % of the time values where the solution will be computed
% if ~all(ismember(p.t_save, t)) 
%     disp("Error: the values where the solution is saved, p.t_save, must be a subset of the values where the solution is computed, t.")
%     status = 1;
%     return
% end

end

