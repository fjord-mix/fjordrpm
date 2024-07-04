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
    p.t_save = t;
end

% Get the indices of the timesteps to save at
inx = find(ismember(t, p.t_save));

s.t = t(inx);

% Box variables
s.H = s.H(:,inx);
s.T = s.T(:,inx);
s.S = s.S(:,inx);
s.V = s.V(:,inx);
s.I = s.I(:,inx);

% Vertical grid
ints = -[0;cumsum(s.H(:,1))];
s.z = 0.5*(ints(1:end-1)+ints(2:end));

% Glacier exchanges
s.QVg = s.QVg(:,inx);
s.QTg = s.QTg(:,inx);
s.QSg = s.QSg(:,inx);

% Shelf exchanges
s.QVs = s.QVs(:,inx);
s.QTs = s.QTs(:,inx);
s.QSs = s.QSs(:,inx);
s.Se = s.Se(:,inx);
s.Te = s.Te(:,inx);
s.phi = s.phi(:,inx);

% Vertical mixing
s.QVk = s.QVk(:,inx);
s.QTk = s.QTk(:,inx);
s.QSk = s.QSk(:,inx);

% Vertical fluxes
s.QVv = s.QVv(:,inx);
s.QTv = s.QTv(:,inx);
s.QSv = s.QSv(:,inx);

% Iceberg fluxes
s.QVi = s.QVi(:,inx);
s.QTi = s.QTi(:,inx);
s.QSi = s.QSi(:,inx);

% Return forcing on same time step as forcings (in results structure to
% prevent overwriting).
s.Ss = f.Ss(:,inx);
s.Ts = f.Ts(:,inx);
s.Qsg = f.Qsg(inx);

s.status = status;

end
