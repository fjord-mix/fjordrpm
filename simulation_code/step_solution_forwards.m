function s = step_solution_forwards(i, p, s)

% STEP_SOLUTION_FOWARDS Does an Euler forward timestep of T and S.
%   s = STEP_SOLUTION_FORWARDS(i, p, s) performs an Euler step for the
%   simulation for parameters structure p and solution s at time step i.

% compute the FjordRPM tracer variables at timestep i+1
% note the plume fluxes are summed over the number of plumes
s.T(:,i+1) = s.T(:,i)+s.dt(i)*p.sid*(sum(s.QTp(:,:,i),1)'+s.QTs(:,i)+s.QTk(:,i)+s.QTi(:,i)+s.QTv(:,i))./s.V;
s.S(:,i+1) = s.S(:,i)+s.dt(i)*p.sid*(sum(s.QSp(:,:,i),1)'+s.QSs(:,i)+s.QSk(:,i)+s.QSi(:,i)+s.QSv(:,i))./s.V;

end