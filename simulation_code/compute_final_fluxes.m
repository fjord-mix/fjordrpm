function s = compute_zmodel_final_fluxes(p, f, s, t)

% COMPUTE_ZMODEL_FINAL_FLUXES Compute and store the fluxes for the solution
% at the final timestep.
%   S =  COMPUTE_ZMODEL_FINAL_FLUXES(P, F, S) computes the final fluxes for
%   input parameters P, forcing structure F andreturns solution structure
%   S.

% Compute the fluxes at the final timestep
[Q, E, I] = compute_zmodel_fluxes(length(t), p, f, s);

% Store the fluxes
% Glacier fluxes
s.QVg(:,end) = Q.QVg;
s.QTg(:,end) = Q.QTg;
s.QSg(:,end) = Q.QSg;

% Shelf fluxes
s.QVs(:,end) = Q.QVs;
s.QTs(:,end) = Q.QTs;
s.QSs(:,end) = Q.QSs;
s.Se(:,end) = E.Se;
s.Te(:,end) = E.Te;
s.phi(:,end) = E.phi;

% Mixing fluxes
s.QVk(:,end) = Q.QVk;
s.QTk(:,end) = Q.QTk;
s.QSk(:,end) = Q.QSk;

% Vertical fluxes
s.QVv(:,end) = Q.QVv;
s.QTv(:,end) = Q.QTv;
s.QSv(:,end) = Q.QSv;

% Iceberg fluxes
s.QIi(:,end) = Q.QIi;
s.QTi(:,end) = Q.QTi;
s.QSi(:,end) = Q.QSi; 
s.QVmi(:,end) = Q.QVmi;
s.M(:,end) = I.M;

end
