function s = compute_zmodel_fluxes(i, p, f, s)


 H =s.H(:,i);
 T = s.T(:,i);
 S = s.S(:,i);
 Qsg = f.Qsg(i);
 zs= f.zs;
 Ts = f.Ts(:,i);
 Ss = f.Ss(:,i);
 V = s.V(:,i);
 I = s.I(:,i);
 zi = f.zi;


% Calculate plume fluxes.
Qg = get_zmodel_plume_fluxes(i, p, f, s);

% Calculate shelf fluxes.
Qs = get_zmodel_shelf_fluxes(i, p, f, s);

% Calculate vertical mixing fluxes.
[QVk ,QTk ,QSk ] = ...
    get_mixing_fluxes(H ,T ,S ,p);

% Calculate iceberg fluxes.
[QIi ,QTi ,QSi ,M , QVmi, QTmi, QSmi] = ...
    get_iceberg_fluxes(H ,T ,S ,I ,p);

% Calculate vertical fluxes.
[QVv ,QTv ,QSv ] = ...
    get_vertical_fluxes(Qg.V-Qs.V+QVmi,T,S,p);

s.QVg(:,i) = Qg.V;
s.QTg(:,i) = Qg.T;
s.QSg(:,i) = Qg.S;

s.QVs(:,i) = Qs.V;
s.QTs(:,i) = Qs.T;
s.QSs(:,i) = Qs.S;
s.Se(:,i) = Qs.Se;
s.Te(:,i) = Qs.Te;
s.phi(:,i) = Qs.phi;

s.QVk(:,i) = QVk;
s.QTk(:,i) = QTk;
s.QSk(:,i) = QSk;
s.QIi(:,i) = QIi;
s.QTi(:,i) = QTi;
s.QSi(:,i) = QSi;
s.M(:,i) = M;
s.QTmi(:,i) = QTmi;
s.Smi(:,i) = QSmi;
s.QVv(:,i) = QVv;
s.QTv(:,i) = QTv;
s.QSv(:,i) = QSv;

end
