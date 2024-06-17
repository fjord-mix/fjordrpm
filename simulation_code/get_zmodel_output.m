function s = get_zmodel_output(p, f, t, s, status)

% Just save daily values as otherwise high time resolution results in large output files
dtdaily = max(1, p.dt);
if status == 0
    int = round(dtdaily/p.dt);
else
    int = 1; % if something went wrong, we want all time steps to properly understand what happened
end
s.t = t(1:int:end-1);

% box variables
s.H = s.H(:,1:int:end-1);
s.T = s.T(:,1:int:end-1);
s.S = s.S(:,1:int:end-1);
s.V = s.V(:,1:int:end-1);
s.I = s.I(:,1:int:end-1);

% glacier exchanges
s.QVg = s.QVg(:,1:int:end);
s.QTg = s.QTg(:,1:int:end);
s.QSg = s.QSg(:,1:int:end);

% shelf exchanges
s.QVs = s.QVs(:,1:int:end);
s.QTs = s.QTs(:,1:int:end);
s.QSs = s.QSs(:,1:int:end);
s.Se = s.Se(:,1:int:end);
s.Te = s.Te(:,1:int:end);
s.phi = s.phi(:,1:int:end);

% vertical mixing
s.QVk = s.QVk(:,1:int:end);
s.QTk = s.QTk(:,1:int:end);
s.QSk = s.QSk(:,1:int:end);

% vertical fluxes
s.QVv = s.QVv(:,1:int:end);
s.QTv = s.QTv(:,1:int:end);
s.QSv = s.QSv(:,1:int:end);

% iceberg fluxes
s.QIi = s.QIi(:,1:int:end);
s.QTi = s.QTi(:,1:int:end);
s.QSi = s.QSi(:,1:int:end);
s.QVmi = s.QVmi(:, 1:int:end);
s.M = s.M(:,1:int:end);

% for iceberg fluxes also calculate and save fjord-integrated values
s.IT = sum(s.I); % fjord iceberg surface area
s.MT = sum(s.QIi); % total iceberg melt flux
% s.ET = p.W*p.L*trapz(f.zi,p.E0*s.I); % total iceberg export flux

% return forcing on same time step as forcings (in results structure to prevent overwriting)
s.Ss = f.Ss(:,1:int:end-1);
s.Ts = f.Ts(:,1:int:end-1);
s.Qsg = f.Qsg(1:int:end-1);
s.D = f.D(1:int:end-1);

s.status = status;
end
