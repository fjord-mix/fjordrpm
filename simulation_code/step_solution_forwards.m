function s = step_solution_forwards(i, p, s)

% STEP_FOWARDS Euler step for the zmodel simulation.
%   S = STEP_FORWARDS(I, P, S) performs an Euler step for the
%   zmodel simulation for parameters structure P and given solution S

% Compute the zmodel tracer variables at timestep i+1 using the timestep
% specified by the input parameters.
s.T(:,i+1) = s.T(:,i)+s.dt(i)*p.sid*(s.QTg(:,i)-s.QTs(:,i)+s.QTk(:,i)+s.QTi(:,i)+s.QTv(:,i))./s.V;
s.S(:,i+1) = s.S(:,i)+s.dt(i)*p.sid*(s.QSg(:,i)-s.QSs(:,i)+s.QSk(:,i)+s.QSi(:,i)+s.QSv(:,i))./s.V;

end
