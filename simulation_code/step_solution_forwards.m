function s = step_solution_forwards(i, p, s)

% STEP_SOLUTION_FOWARDS Does an Euler forward timestep of T and S.
%   s = STEP_SOLUTION_FORWARDS(i, p, s) performs an Euler step for the
%   simulation for parameters structure p and solution s at time step i.

% compute the FjordRPM tracer variables at timestep i+1
% case with multiple plumes needs special treatment
if size(s.QTp,1)==1
    s.T(:,i+1) = s.T(:,i)+s.dt(i)*p.sid*(s.QTp(:,:,i)'+s.QTs(:,i)+s.QTk(:,i)+s.QTi(:,i)+s.QTv(:,i))./s.V;
    s.S(:,i+1) = s.S(:,i)+s.dt(i)*p.sid*(s.QSp(:,:,i)'+s.QSs(:,i)+s.QSk(:,i)+s.QSi(:,i)+s.QSv(:,i))./s.V;
% or if only 1 plume
else
    s.T(:,i+1) = s.T(:,i)+s.dt(i)*p.sid*(sum(s.QTp(:,:,i))'+s.QTs(:,i)+s.QTk(:,i)+s.QTi(:,i)+s.QTv(:,i))./s.V;
    s.S(:,i+1) = s.S(:,i)+s.dt(i)*p.sid*(sum(s.QSp(:,:,i))'+s.QSs(:,i)+s.QSk(:,i)+s.QSi(:,i)+s.QSv(:,i))./s.V;
end

end