function s = get_final_output(p, f, t, s, status)

% GET_OUTPUT Get the output of the z-model simulation.
%   S = GET_OUTPUT(P, F, T, S, STATUS) gets the zmodel output for input
%   parameters P, forcing structure F, time vector T and solution structure
%   S, and returns solution structure S on timestepping specified by the
%   user, including forcing in the same time steps as S.

% Save values to output file as specified by the input parameter, unless
% something went wrong, then we want all time steps to properly understand
% what happened.
if status == 1
    p.t_save = t+1; % Indexing must start from 1 not 0
end

s.t = t(p.t_save+1);

% Box variables
s.H = s.H(:,p.t_save+1);
s.T = s.T(:,p.t_save+1);
s.S = s.S(:,p.t_save+1);
s.V = s.V(:,p.t_save+1);
s.I = s.I(:,p.t_save+1);

% Glacier exchanges
s.QVg = s.QVg(:,p.t_save+1);
s.QTg = s.QTg(:,p.t_save+1);
s.QSg = s.QSg(:,p.t_save+1);

% Shelf exchanges
s.QVs = s.QVs(:,p.t_save+1);
s.QTs = s.QTs(:,p.t_save+1);
s.QSs = s.QSs(:,p.t_save+1);
s.Se = s.Se(:,p.t_save+1);
s.Te = s.Te(:,p.t_save+1);
s.phi = s.phi(:,p.t_save+1);

% Vertical mixing
s.QVk = s.QVk(:,p.t_save+1);
s.QTk = s.QTk(:,p.t_save+1);
s.QSk = s.QSk(:,p.t_save+1);

% Vertical fluxes
s.QVv = s.QVv(:,p.t_save+1);
s.QTv = s.QTv(:,p.t_save+1);
s.QSv = s.QSv(:,p.t_save+1);

% Iceberg fluxes
s.QIi = s.QIi(:,p.t_save+1);
s.QTi = s.QTi(:,p.t_save+1);
s.QSi = s.QSi(:,p.t_save+1); 
s.QVmi = s.QVmi(:,p.t_save+1);
s.M = s.M(:,p.t_save+1);

% For iceberg fluxes, also calculate and save fjord-integrated values.
s.IT = sum(s.I); % fjord iceberg surface area
s.MT = sum(s.QIi); % total iceberg melt flux
% s.ET = p.W*p.L*trapz(f.zi,p.E0*s.I); % total iceberg export flux

% Return forcing on same time step as forcings (in results structure to
% prevent overwriting).
s.Ss = f.Ss(:,p.t_save+1);
s.Ts = f.Ts(:,p.t_save+1);
s.Qsg = f.Qsg(p.t_save+1);
s.D = f.D(p.t_save+1);

s.status = status;
end
