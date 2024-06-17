function Q = compute_zmodel_fluxes(i, p, f, s)

% Calculate plume fluxes.
[Q.QVg, Q.QTg, Q.QSg] = get_zmodel_plume_fluxes(i, p, f, s);

% Calculate shelf fluxes.
[Q.QVs, Q.QTs, Q.QSs, Q.Se, Q.Te, Q.phi ] = get_zmodel_shelf_fluxes(i, p, f, s);

% Calculate vertical mixing fluxes.
[Q.QVk, Q.QTk, Q.QSk] = get_zmodel_mixing_fluxes(i, p, s);

% Calculate iceberg fluxes.
[Q.QIi, Q.QTi, Q.QSi, Q.M, Q.QVmi, Q.QTmi, Q.QSmi] = get_zmodel_iceberg_fluxes(i, p, s);

% Calculate vertical fluxes.
[Q.QVv, Q.QTv, Q.QSv ] = get_zmodel_vertical_fluxes(i, p, s, Q.QVg-Q.QVs+Q.QVmi );

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


end
