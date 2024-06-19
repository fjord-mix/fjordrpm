function Tr = step_zmodel_forwards(i, p, s, Q)

% STEP_ZMODEL_FOWARDS Euler step for the zmodel simulation.
%   TR = STEP_ZMODEL_FORWARDS(I, P, S, Q) performs an Euler step for the
%   zmodel simulation for parameters structure P and given solution S and
%   returns updated tracer variables stored in a structure TR.

% Compute the zmodel tracer variables at timestep i+1.
Tr.V = s.V(:,i) +p.dt*p.sid*(Q.QVg-Q.QVs+Q.QVk+Q.QVmi+Q.QVv);
Tr.VT = s.VT(:,i)+p.dt*p.sid*(Q.QTg-Q.QTs+Q.QTk+Q.QTi+Q.QTmi+Q.QTv);
Tr.VS = s.VS(:,i)+p.dt*p.sid*(Q.QSg-Q.QSs+Q.QSk+Q.QSi+Q.QSmi+Q.QSv);
Tr.H = Tr.V/(p.W*p.L);
Tr.T = Tr.VT./Tr.V;
Tr.S = Tr.VS./Tr.V;

% No update to static icebergs.
Tr.I = s.I(:,i);
%     I = I+
%     (1-p.icestatic)*p.dt*p.sid*((D/(p.W*p.L))*f.xi-M.*I-p.uIce/p.L*I);

end
