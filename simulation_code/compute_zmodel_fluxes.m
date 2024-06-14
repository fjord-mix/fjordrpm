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
Qg_i = get_zmodel_plume_fluxes(i, p, f, s);

% Calculate shelf fluxes.
[QVs ,QTs ,QSs ,Se ,Te ,phi ] = ...
    get_shelf_fluxes(H ,T ,S ,zs,Ts ,Ss ,Qsg,p, s);

% Calculate vertical mixing fluxes.
[QVk ,QTk ,QSk ] = ...
    get_mixing_fluxes(H ,T ,S ,p);

% % Calculate "artificial" fluxes.
% [QVb ,QTb ,QSb ] = ...
%     get_artificial_fluxes(QVg -QVs +QVk ,H ,V ,T ,S ,zs,Ts ,Ss ,p);

% Calculate iceberg fluxes.
[QIi ,QTi ,QSi ,M , QVmi, QTmi, QSmi] = ...
    get_iceberg_fluxes(H ,T ,S ,I ,p);

% Calculate vertical fluxes.
[QVv ,QTv ,QSv ] = ...
    get_vertical_fluxes(Qg_i.V-QVs+QVmi,T,S,p);

s.QVg(:,i) = Qg_i.V;
s.QTg(:,i) = Qg_i.T;
s.QSg(:,i) = Qg_i.S;

s.QVs(:,i) = QVs;
s.QTs(:,i) = QTs;
s.QSs(:,i) = QSs;
s.Se(:,i) = Se;
s.Te(:,i) = Te;
s.phi(:,i) = phi;
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
