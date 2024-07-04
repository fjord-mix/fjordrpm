function Tr = step_solution_forwards(i, p, s, Q)

% STEP_FOWARDS Euler step for the zmodel simulation.
%   TR = STEP_FORWARDS(I, P, S, Q) performs an Euler step for the
%   zmodel simulation for parameters structure P and given solution S and
%   returns updated tracer variables stored in a structure TR.

% Compute the zmodel tracer variables at timestep i+1 using the timestep
% specified by the input parameters.
Tr.T = s.T(:,i)+s.dt(i)*p.sid*(Q.QTg-Q.QTs+Q.QTk+Q.QTi+Q.QTv)./s.V;
Tr.S = s.S(:,i)+s.dt(i)*p.sid*(Q.QSg-Q.QSs+Q.QSk+Q.QSi+Q.QSv)./s.V;

end
