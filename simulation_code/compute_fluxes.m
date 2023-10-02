function [QVg, QTg, QSg, QVs,QTs,QSs,Se,Te,phi ,QVk ,QTk ,QSk , QVb ,...
    QTb ,QSb , QIi ,QTi ,QSi ,M  ]...
    = compute_fluxes (H ,T ,S ,Qsg, p, zs, Ts ,Ss, V , I ,zi)

% Calculate plume fluxes.
[QVg ,QTg ,QSg ] = ...
    get_plume_fluxes(H ,T ,S ,Qsg,p);

% Calculate shelf fluxes.
[QVs ,QTs ,QSs ,Se ,Te ,phi ] = ...
    get_shelf_fluxes(H ,T ,S ,zs,Ts ,Ss ,Qsg,p);

% Calculate vertical mixing fluxes.
[QVk ,QTk ,QSk ] = ...
    get_mixing_fluxes(H ,T ,S ,QVg ,QVs ,p);

% Calculate "artificial" fluxes.
[QVb ,QTb ,QSb ] = ...
    get_artificial_fluxes(QVg -QVs +QVk ,H ,V ,T ,S ,zs,Ts ,Ss ,p);

% Calculate iceberg fluxes.
[QIi ,QTi ,QSi ,M ] = ...
    get_iceberg_fluxes(H ,T ,S ,I ,zi,p);
end