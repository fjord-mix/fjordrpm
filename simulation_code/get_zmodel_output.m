function s = get_zmodel_output(p, f, t, H, T, S, V, I, ...
    QVg, QTg, QSg, QVs, QTs, QSs, Se, Te, phi, QVk, QTk, QSk, ...
    QVv, QTv, QSv, QIi, QTi, QSi, QVmi, M)

% Just save ~daily values
% as otherwise high time resolution results in large output files
dtdaily = max(1, p.dt);
if s.status == 0
    int = round(dtdaily/p.dt);
else
    int = 1; % if something went wrong, we want all time steps to properly understand what happened
end
int = round(dtdaily/p.dt);
s.t = t(1:int:end-1);

% box variables
s.H = H(:,1:int:end-1);
s.T = T(:,1:int:end-1);
s.S = S(:,1:int:end-1);
s.V = V(:,1:int:end-1);
s.I = I(:,1:int:end-1);

% glacier exchanges
s.QVg = QVg(:,1:int:end);
s.QTg = QTg(:,1:int:end);
s.QSg = QSg(:,1:int:end);

% shelf exchanges
s.QVs = QVs(:,1:int:end);
s.QTs = QTs(:,1:int:end);
s.QSs = QSs(:,1:int:end);
s.Se = Se(:,1:int:end);
s.Te = Te(:,1:int:end);
s.phi = phi(:,1:int:end);

% vertical mixing
s.QVk = QVk(:,1:int:end);
s.QTk = QTk(:,1:int:end);
s.QSk = QSk(:,1:int:end);

% vertical fluxes
s.QVv = QVv(:,1:int:end);
s.QTv = QTv(:,1:int:end);
s.QSv = QSv(:,1:int:end);

% iceberg fluxes
s.QIi = QIi(:,1:int:end);
s.QTi = QTi(:,1:int:end);
s.QSi = QSi(:,1:int:end);
s.QVmi = QVmi(:, 1:int:end);
s.M = M(:,1:int:end);

% for iceberg fluxes also calculate and save fjord-integrated values
s.IT = p.W*p.L*trapz(f.zi,s.I); % fjord iceberg volume
s.MT = p.W*p.L*trapz(f.zi,s.M(:,1:size(s.I,2)).*s.I); % total iceberg melt flux
s.ET = p.W*p.L*trapz(f.zi,p.E0*s.I); % total iceberg export flux

% return forcing on same time step as forcings (in results structure to prevent overwriting)
s.Ss = f.Ss(:,1:int:end-1);
s.Ts = f.Ts(:,1:int:end-1);
s.Qsg = f.Qsg(1:int:end-1);
s.D = f.D(1:int:end-1);

end
