function s = store_zmodel_solution(i, s, Q_i, Z_ip1)

% Store zmodel fluxes at timestep i
s.QVg(:,i) = Q_i.QVg;
s.QTg(:,i) = Q_i.QTg;
s.QSg(:,i) = Q_i.QSg;

s.QVs(:,i) = Q_i.QVs;
s.QTs(:,i) = Q_i.QTs;
s.QSs(:,i) = Q_i.QSs;
s.Se(:,i) = Q_i.Se;
s.Te(:,i) = Q_i.Te;
s.phi(:,i) =Q_i.phi;

s.QVk(:,i) = Q_i.QVk;
s.QTk(:,i) = Q_i.QTk;
s.QSk(:,i) = Q_i.QSk;

s.QIi(:,i) = Q_i.QIi;
s.QTi(:,i) = Q_i.QTi;
s.QSi(:,i) = Q_i.QSi;
s.M(:,i) = Q_i.M;
s.QVmi(:,i) = Q_i.QVmi;
s.QTmi(:,i) = Q_i.QTmi;
s.Smi(:,i) = Q_i.QSmi;

s.QVv(:,i) = Q_i.QVv;
s.QTv(:,i) = Q_i.QTv;
s.QSv(:,i) = Q_i.QSv;

% Store zmodel variables at timestep i+1
s.V(:,i+1) = Z_ip1.V;
s.T(:,i+1) = Z_ip1.T;
s.S(:,i+1) = Z_ip1.S;
s.H(:,i+1) = Z_ip1.H;
s.VT(:,i+1) = Z_ip1.VT;
s.VS(:,i+1) = Z_ip1.VS;

end
