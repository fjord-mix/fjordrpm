function s = store_zmodel_solution(i, s, Q, E, I, Tr)

% STORE_ZMODEL_SOLUTION z-model simulation.
%   S = STORE_ZMODEL_SOLUTION(I, S, Q, E, I, TR) stores the zmodel solution
%   Q, E, I, TR after iteration I in structure S.

% Store zmodel fluxes and shelf/ice variables at timestep i.
s.QVg(:,i) = Q.QVg;
s.QTg(:,i) = Q.QTg;
s.QSg(:,i) = Q.QSg;

s.QVs(:,i) = Q.QVs;
s.QTs(:,i) = Q.QTs;
s.QSs(:,i) = Q.QSs;
s.Se(:,i) = E.Se;
s.Te(:,i) = E.Te;
s.phi(:,i) =E.phi;

s.QVk(:,i) = Q.QVk;
s.QTk(:,i) = Q.QTk;
s.QSk(:,i) = Q.QSk;

s.QIi(:,i) = Q.QIi;
s.QTi(:,i) = Q.QTi;
s.QSi(:,i) = Q.QSi;
s.QVmi(:,i) = Q.QVmi;
s.QTmi(:,i) = Q.QTmi;
s.Smi(:,i) = Q.QSmi;
s.M(:,i) = I.M;

s.QVv(:,i) = Q.QVv;
s.QTv(:,i) = Q.QTv;
s.QSv(:,i) = Q.QSv;

% Store zmodel tracer variables at timestep i+1.
s.V(:,i+1) = Tr.V;
s.T(:,i+1) = Tr.T;
s.S(:,i+1) = Tr.S;
s.H(:,i+1) = Tr.H;
s.VT(:,i+1) = Tr.VT;
s.VS(:,i+1) = Tr.VS;

end
