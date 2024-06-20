function s = get_zmodel_output(p, f, t, s, status)

% GET_ZMODEL_OUTPUT Get the output of the z-model simulation.
%   S = GET_ZMODEL_OUTPUT(P, F, T, S, STATUS) gets the zmodel output for
%   input parameters P, forcing structure F, time vector T and solution
%   structure S, and returns solution structure S on sparser timestepping
%   if the input timestep is less than one day, including forcing in the
%   same time steps as S.

% Save daily or sparser values, as otherwise high time resolution results
% in large output files.
dtdaily = max(1, s.dt);
if status == 0
    int = round(dtdaily/s.dt);
else
    % If something went wrong, we want all time steps to properly
    % understand what happened.
    int = s.dt;
end
s.t = t(1:int:end);

% Box variables
s.H = s.H(:,1:int:end);
s.T = s.T(:,1:int:end);
s.S = s.S(:,1:int:end);
s.V = s.V(:,1:int:end);
s.I = s.I(:,1:int:end);

% Glacier exchanges
s.QVg = s.QVg(:,1:int:end);
s.QTg = s.QTg(:,1:int:end);
s.QSg = s.QSg(:,1:int:end);

% Shelf exchanges
s.QVs = s.QVs(:,1:int:end);
s.QTs = s.QTs(:,1:int:end);
s.QSs = s.QSs(:,1:int:end);
s.Se = s.Se(:,1:int:end);
s.Te = s.Te(:,1:int:end);
s.phi = s.phi(:,1:int:end);

% Vertical mixing
s.QVk = s.QVk(:,1:int:end);
s.QTk = s.QTk(:,1:int:end);
s.QSk = s.QSk(:,1:int:end);

% Vertical fluxes
s.QVv = s.QVv(:,1:int:end);
s.QTv = s.QTv(:,1:int:end);
s.QSv = s.QSv(:,1:int:end);

% Iceberg fluxes
s.QIi = s.QIi(:,1:int:end);
s.QTi = s.QTi(:,1:int:end);
s.QSi = s.QSi(:,1:int:end); 
s.QVmi = s.QVmi(:,1:int:end);
s.M = s.M(:,1:int:end);

% For iceberg fluxes, also calculate and save fjord-integrated values.
s.IT = sum(s.I); % fjord iceberg surface area
s.MT = sum(s.QIi); % total iceberg melt flux
% s.ET = p.W*p.L*trapz(f.zi,p.E0*s.I); % total iceberg export flux

% Return forcing on same time step as forcings (in results structure to
% prevent overwriting).
s.Ss = f.Ss(:,1:int:end);
s.Ts = f.Ts(:,1:int:end);
s.Qsg = f.Qsg(1:int:end);
s.D = f.D(1:int:end);

s.status = status;
end
