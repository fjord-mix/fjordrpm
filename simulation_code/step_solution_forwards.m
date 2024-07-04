function Tr = step_solution_forwards(i, p, s, Q)

% STEP_FOWARDS Euler step for the zmodel simulation.
%   TR = STEP_FORWARDS(I, P, S, Q) performs an Euler step for the
%   zmodel simulation for parameters structure P and given solution S and
%   returns updated tracer variables stored in a structure TR.

% No update to static volume, thickness, or icebergs.
Tr.I = s.I(:,i);
Tr.V = s.V(:,i);
Tr.H = Tr.V/(p.W*p.L);

% Compute the zmodel tracer variables at timestep i+1 using the timestep
% specified by the input parameters.
Tr.VT = s.VT(:,i)+s.dt(i)*p.sid*(Q.QTg-Q.QTs+Q.QTk+Q.QTi+Q.QTv);
Tr.VS = s.VS(:,i)+s.dt(i)*p.sid*(Q.QSg-Q.QSs+Q.QSk+Q.QSi+Q.QSv);
Tr.T = Tr.VT./Tr.V;
Tr.S = Tr.VS./Tr.V;

end
