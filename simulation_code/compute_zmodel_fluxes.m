function s = compute_zmodel_fluxes(i, p, f, s)

% Calculate plume fluxes.
Qg = get_zmodel_plume_fluxes(i, p, f, s);

% Calculate shelf fluxes.
Qs = get_zmodel_shelf_fluxes(i, p, f, s);

% Calculate vertical mixing fluxes.
Qk = get_zmodel_mixing_fluxes(i, p, s);

% Calculate iceberg fluxes.
Qi = get_zmodel_iceberg_fluxes(i, p, s);

% Calculate vertical fluxes.
Qv = get_zmodel_vertical_fluxes(i, p, s, Qg, Qs, Qi);

s.QVg(:,i) = Qg.V;
s.QTg(:,i) = Qg.T;
s.QSg(:,i) = Qg.S;

s.QVs(:,i) = Qs.V;
s.QTs(:,i) = Qs.T;
s.QSs(:,i) = Qs.S;
s.Se(:,i) = Qs.Se;
s.Te(:,i) = Qs.Te;
s.phi(:,i) = Qs.phi;

s.QVk(:,i) = Qk.V;
s.QTk(:,i) = Qk.T;
s.QSk(:,i) = Qk.S;

s.QIi(:,i) = Qi.I;
s.QTi(:,i) = Qi.T;
s.QSi(:,i) = Qi.S;
s.M(:,i) = Qi.M;
s.QVmi(:,i) = Qi.Vm;
s.QTmi(:,i) = Qi.Tm;
s.Smi(:,i) = Qi.Sm;

s.QVv(:,i) = Qv.V;
s.QTv(:,i) = Qv.T;
s.QSv(:,i) = Qv.S;

end
