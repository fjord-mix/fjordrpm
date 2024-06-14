function [Q, s] = compute_zmodel_fluxes(i, p, f, s)

% Calculate plume fluxes.
Q.Qg = get_zmodel_plume_fluxes(i, p, f, s);

% Calculate shelf fluxes.
Q.Qs = get_zmodel_shelf_fluxes(i, p, f, s);

% Calculate vertical mixing fluxes.
Q.Qk = get_zmodel_mixing_fluxes(i, p, s);

% Calculate iceberg fluxes.
Q.Qi = get_zmodel_iceberg_fluxes(i, p, s);

% Calculate vertical fluxes.
Q.Qv = get_zmodel_vertical_fluxes(i, p, s, Q);

s.QVg(:,i) = Q.Qg.V;
s.QTg(:,i) = Q.Qg.T;
s.QSg(:,i) = Q.Qg.S;

s.QVs(:,i) = Q.Qs.V;
s.QTs(:,i) = Q.Qs.T;
s.QSs(:,i) = Q.Qs.S;
s.Se(:,i) = Q.Qs.Se;
s.Te(:,i) = Q.Qs.Te;
s.phi(:,i) =Q.Qs.phi;

s.QVk(:,i) = Q.Qk.V;
s.QTk(:,i) = Q.Qk.T;
s.QSk(:,i) = Q.Qk.S;

s.QIi(:,i) = Q.Qi.I;
s.QTi(:,i) = Q.Qi.T;
s.QSi(:,i) = Q.Qi.S;
s.M(:,i) = Q.Qi.M;
s.QVmi(:,i) = Q.Qi.Vm;
s.QTmi(:,i) = Q.Qi.Tm;
s.Smi(:,i) = Q.Qi.Sm;

s.QVv(:,i) = Q.Qv.V;
s.QTv(:,i) = Q.Qv.T;
s.QSv(:,i) = Q.Qv.S;

end
