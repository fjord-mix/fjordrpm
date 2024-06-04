function [V, T, S, I, VT, VS,...
    QVg, QTg, QSg, QVs, QTs, QSs,...
    QVk, QTk, QSk, QVmi, QTmi, QSmi,QIi,QTi,QSi,...
    QVv,QTv,QSv,Te,Se, phi, M] = initialise_zmodel_variables(p, f, a, t)

[T, S] = deal(zeros(p.N, length(t)));

[QVg,QTg,QSg,...
    QVs,QTs,QSs,...
    QVk,QTk,QSk,...
    QVmi,QTmi,QSmi,...
    QIi,QTi,QSi,...
    QVv,QTv,QSv,...
    Te,Se, phi] = deal(zeros(p.N,length(t)-1));

[I, M] = deal(zeros(length(f.zi),length(t)-1));

%% Initialise variables according with the boundary and initial conditions.

V = a.H0'*p.W*p.L; % volume of layers
T(:,1) = a.T0; % temperature
S(:,1) = a.S0; % salinity
I(:,1) = a.I0; % iceberg concentration
VT(:,1) = V.*T(:,1); % heat content
VS(:,1) = V.*S(:,1); % salt content

end