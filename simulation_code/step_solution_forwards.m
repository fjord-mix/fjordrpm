function s = step_solution_forwards(i, p, s)

% STEP_SOLUTION_FOWARDS Euler timestep for the simulation.
%   s = STEP_SOLUTION_FORWARDS(i, p, s) performs an Euler step for the
%   simulation for parameters structure p and given solution s.

% Compute the zmodel tracer variables at timestep i+1
s.T(:,i+1) = s.T(:,i)+s.dt(i)*p.sid*(s.QTg(:,i)+s.QTs(:,i)+s.QTk(:,i)+s.QTi(:,i)+s.QTv(:,i))./s.V;
s.S(:,i+1) = s.S(:,i)+s.dt(i)*p.sid*(s.QSg(:,i)+s.QSs(:,i)+s.QSk(:,i)+s.QSi(:,i)+s.QSv(:,i))./s.V;

end
