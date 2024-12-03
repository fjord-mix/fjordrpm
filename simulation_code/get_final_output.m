function s = get_final_output(p, t, s, status)

% GET_FINAL_OUTPUT Get the output of the simulation.
%   s = GET_FINAL_OUTPUT(p, t, s, status) gets the model output for
%   input parameters p, time vector t and solution structure s, and 
%   returns solution structure s on timestepping specified by the user.

% Compute fluxes one last time to get same number of values as for tracers
s = compute_fluxes(size(s.T,2), p, s);

% Save values to output file as specified by the input parameter, unless
% (a) something went wrong, because we want all time steps to properly understand
% what happened, or (b) we did not specify particular time steps for saving
if status == 1 || ~isfield(p,'t_save')
    p.t_save = t;
end

% Get the indices of the timesteps to save at
inx = find(ismember(t, p.t_save));

s.t = t(inx);

% Layer variables
s.T = s.T(:,inx);
s.S = s.S(:,inx);

% Vertical grid
ints = -[0;cumsum(s.H(:,1))];
s.z = 0.5*(ints(1:end-1)+ints(2:end));

% Plume exchange
s.QVp = s.QVp(:,:,inx);
s.QTp = s.QTp(:,:,inx);
s.QSp = s.QSp(:,:,inx);
s.QMp = s.QMp(:,:,inx);
for j=1:length(p.wp)
    s.plumemeltrate(j,:,:) = p.sid*squeeze(s.QMp(j,:,:))./(p.wp(j)*s.H);
end
s.knb = s.knb(inx);

% Shelf exchange
s.QVs = s.QVs(:,inx);
s.QTs = s.QTs(:,inx);
s.QSs = s.QSs(:,inx);
s.Ss = s.Ss(:,inx);
s.Ts = s.Ts(:,inx);
s.phi = s.phi(:,inx);

% Vertical mixing
s.QVk = s.QVk(:,inx);
s.QTk = s.QTk(:,inx);
s.QSk = s.QSk(:,inx);

% Vertical fluxes
s.QVv = s.QVv(:,inx);
s.QTv = s.QTv(:,inx);
s.QSv = s.QSv(:,inx);

% Iceberg fluxes and melt rates
s.QVi = s.QVi(:,inx);
s.QTi = s.QTi(:,inx);
s.QSi = s.QSi(:,inx);
s.QMi = s.QMi(:,inx);
s.icebergmeltrate = p.sid*s.QMi./s.I;
s.icebergmeltrate(s.I==0,:) = 0;

% Subglacial discharge
s.Qsg = s.Qsg(:,inx);


s.status = status;

end
