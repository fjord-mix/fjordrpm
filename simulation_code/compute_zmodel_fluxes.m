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
[Q.QVv, Q.QTv, Q.QSv ] = get_zmodel_vertical_fluxes(i, p, s, Q.QVg-Q.QVs+Q.QVmi);

end
