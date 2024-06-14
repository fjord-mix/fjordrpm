function [s, Q] = compute_zmodel_fluxes(i, p, f, s)

% Calculate plume fluxes.
[Q.QVg ,Q.QTg ,Q.QSg] = get_zmodel_plume_fluxes(i, p, f, s);

% Calculate shelf fluxes.
[Q.QVs ,Q.QTs ,Q.QSs ,Q.Se ,Q.Te ,Q.phi ] = get_zmodel_shelf_fluxes(i, p, f, s);

% Calculate vertical mixing fluxes.
[Q.QVk ,Q.QTk ,Q.QSk ] = get_zmodel_mixing_fluxes(i, p, s);

% Calculate iceberg fluxes.
[Q.QIi ,Q.QTi ,Q.QSi ,Q.M , Q.QVmi, Q.QTmi, Q.QSmi] = get_zmodel_iceberg_fluxes(i, p, s);

% Calculate vertical fluxes.
[Q.QVv ,Q.QTv ,Q.QSv ] = get_zmodel_vertical_fluxes(i, p, s, Q.QVg-Q.QVs+Q.QVmi );

% Store zmodel fluxes from this timestep
s.QVg(:,i) = Q.QVg;
s.QTg(:,i) = Q.QTg;
s.QSg(:,i) = Q.QSg;

s.QVs(:,i) = Q.QVs;
s.QTs(:,i) = Q.QTs;
s.QSs(:,i) = Q.QSs;
s.Se(:,i) = Q.Se;
s.Te(:,i) = Q.Te;
s.phi(:,i) =Q.phi;

s.QVk(:,i) = Q.QVk;
s.QTk(:,i) = Q.QTk;
s.QSk(:,i) = Q.QSk;

s.QIi(:,i) = Q.QIi;
s.QTi(:,i) = Q.QTi;
s.QSi(:,i) = Q.QSi;
s.M(:,i) = Q.M;
s.QVmi(:,i) = Q.QVmi;
s.QTmi(:,i) = Q.QTmi;
s.Smi(:,i) = Q.QSmi;

s.QVv(:,i) = Q.QVv;
s.QTv(:,i) = Q.QTv;
s.QSv(:,i) = Q.QSv;

end
