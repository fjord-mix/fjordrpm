function [s,f] = boxmodel(p,f,a,t)

% BOXMODEL Box model simulation.
%   [S,F] = BOXMODEL(P,F,A,T) runs the box model simulation for constant parameters structure P,
%   forcings structure F, initial conditions structure I and time variable
%   T and returns solution structure S and forcing structure F.

s.status = 0; % initial assumption that the model ran successfully

%% Initialise variables
H(:,1) = a.H0; % thickness
T(:,1) = a.T0; % temperature
S(:,1) = a.S0; % salinity
I(:,1) = a.I0; % iceberg concentration
V(:,1) = H(:,1).*double(p.W).*double(p.L); % volume
VT(:,1) = V(:,1).*T(:,1); % heat
VS(:,1) = V(:,1).*S(:,1); % salt

% Preallocate variables
[QVg,QTg,QSg,QVs,QTs,QSs,QVk,QTk,QSk,QVb,QTb,QSb,QIi,QTi,QSi,Te,Se] = deal(zeros(size(H,1),length(t)-1));
phi = zeros(size(H,1)-p.sill,length(t)-1);
M = zeros(length(f.zi),length(t)-1);

%% Error checks
% Check shelf oscillation parameters have been set up correctly.
if (isfield(p,'zd') && p.zd > 0) && (isfield(p,'tw') && p.tw  <= 0)
    disp('Error: must have positive oscillation period if oscillation strength is set.')
    s.status = 1;
    return
end

% Check initialisation is consistent with specified number of layers.
if any([length(H) ~= p.N+p.sill,length(T) ~= p.N+p.sill,length(S) ~= p.N+p.sill])
    disp('Error: Initial conditions not consistent with number of layers');
    s.status = 1; % status == 1 means there was an error
    return
end

% Check bottom box is consistent with sill depth.
if p.sill == 1 && H(end,1) ~= p.H-abs(p.silldepth)
    disp('Error: when p.sill=1, bottom box must have thickness p.H-p.silldepth');
    s.status = 1; % status == 1 means there was an error
    return
end

% Check sum of layer thicknesses is equal to fjord depth.
if abs(sum(H(:,1))-p.H) > 1e-10
    disp('Error: box thicknesses must sum to fjord depth');
    s.status = 1; % status == 1 means there was an error
    return
end

% If layer nudging active, check we have the required nudging inputs.
if ~isnan(p.trelax) && length(p.Snudge) < p.N-1
    disp('Error: incorrect number of nudging values');
    s.status = 1; % status == 1 means there was an error
    return
end

%% Increase the profile resolution for better accuracy when integrating
nz_orig = 1:length(f.zs);
nz_hr = linspace(1,length(f.zs),100*length(f.zs));
zs_hr = interp1(nz_orig,f.zs,nz_hr);
Ts_hr = interp1(f.zs,f.Ts,zs_hr);
Ss_hr = interp1(f.zs,f.Ss,zs_hr);
f.zs = zs_hr;
f.Ts = Ts_hr;
f.Ss = Ss_hr;

%%
if p.plot_runtime
    % hf_track = monitor_boxmodel([],1,H,T,S,f);
    % hf_track = show_boxmodel([],1,t,H,T,S,[],[],[],[],f);
    s_bnds = [min(f.Ss(:)) max(f.Ss(:))+0.1];
end

%% The main loop
for i = 1:length(t)-1

    % Calculate plume fluxes.
    [QVg(:,i),QTg(:,i),QSg(:,i)] = ...
        get_plume_fluxes(H(:,i),T(:,i),S(:,i),f.Qsg(i),p);

    % Calculate shelf fluxes.
    [QVs(:,i),QTs(:,i),QSs(:,i),Se(:,i),Te(:,i),phi(:,i)] = ...
        get_shelf_fluxes(H(:,i),T(:,i),S(:,i),f.zs,f.Ts(:,i),f.Ss(:,i),f.Qsg(i),p);

    % Calculate vertical mixing fluxes.
    [QVk(:,i),QTk(:,i),QSk(:,i)] = ...
        get_mixing_fluxes(H(:,i),T(:,i),S(:,i),QVg(:,i),QVs(:,i),p);

    % Calculate "artificial" fluxes.
    [QVb(:,i),QTb(:,i),QSb(:,i)] = ...
        get_artificial_fluxes(QVg(:,i)-QVs(:,i)+QVk(:,i),H(:,i),V(:,i),T(:,i),S(:,i),f.zs,f.Ts(:,i),f.Ss(:,i),p);

    % Calculate iceberg fluxes.
    [QIi(:,i),QTi(:,i),QSi(:,i),M(:,i)] = ...
        get_iceberg_fluxes(H(:,i),T(:,i),S(:,i),I(:,i),f.zi,p);

    % Step fjord forwards.
    dt = t(i+1)-t(i);
    V(:,i+1)  = V(:,i)+dt*p.sid*(QVg(:,i)-QVs(:,i)+QVk(:,i)+QVb(:,i));
    VT(:,i+1) = VT(:,i)+dt*p.sid*(QTg(:,i)-QTs(:,i)+QTk(:,i)+QTb(:,i)+QTi(:,i));
    VS(:,i+1) = VS(:,i)+dt*p.sid*(QSg(:,i)-QSs(:,i)+QSk(:,i)+QSb(:,i)+QSi(:,i));

    % Calculate thicknesses and tracers.
    H(:,i+1) = V(:,i+1)/(p.W*p.L);
    T(:,i+1) = VT(:,i+1)./V(:,i+1);
    S(:,i+1) = VS(:,i+1)./V(:,i+1);

    % Step icebergs forwards.
    I(:,i+1) = I(:,i)+dt*p.sid*((f.D(i)/(p.W*p.L))*f.xi-M(:,i).*I(:,i)-p.E0*I(:,i));

    % real-time nudging, i.e., nudging values are updated to mimic the current shelf conditions
    if ~isnan(p.trelax) && p.real_time_nudge
        p.Snudge = get_interface_salinities(f.zs,f.Ts(:,i),f.Ss(:,i),p);
    end

    % Plot model evolution (mainly debugging).
    if p.plot_runtime
        % hf_track = monitor_boxmodel(hf_track,i,H,T,S,f);
        % hf_track = show_box_model(hf_track,i,t,H,T,S,QVs,QVg,QVk,QVb,f);
        plot_debug_profile(i,t,f,p,H,S,s_bnds);
    end

    % Break from loop if any layer thickness goes below p.Hmin.
    if ~isempty(find(H(:,i+1) < p.Hmin,1))
        disp('Error: layer thickness dropped below p.Hmin');
        s.status = 1; % status == 1 means there was an error
        break
    end

end

%% Put solution into output structure
% but just save ~daily values
% as otherwise high time resolution results in large output files
dtdaily = 1;
if s.status == 0
    int = round(dtdaily/dt);
else
    int = 1; % if something went wrong, we want all time steps to properly understand what happened
end
s.t = t(1:int:end-1);

% box variables
s.H = H(:,1:int:end-1);
s.T = T(:,1:int:end-1);
s.S = S(:,1:int:end-1);
s.V = V(:,1:int:end-1);
s.I = I(:,1:int:end-1);

% glacier exchanges
s.QVg = QVg(:,1:int:end);
s.QTg = QTg(:,1:int:end);
s.QSg = QSg(:,1:int:end);

% shelf exchanges
s.QVs = QVs(:,1:int:end);
s.QTs = QTs(:,1:int:end);
s.QSs = QSs(:,1:int:end);
s.Se = Se(:,1:int:end);
s.Te = Te(:,1:int:end);
s.phi = phi(:,1:int:end);

% vertical mixing
s.QVk = QVk(:,1:int:end);
s.QTk = QTk(:,1:int:end);
s.QSk = QSk(:,1:int:end);

% artificial fluxes
s.QVb = QVb(:,1:int:end);
s.QTb = QTb(:,1:int:end);
s.QSb = QSb(:,1:int:end);

% iceberg fluxes
s.QIi = QIi(:,1:int:end);
s.QTi = QTi(:,1:int:end);
s.QSi = QSi(:,1:int:end);
s.M = M(:,1:int:end);

% for iceberg fluxes also calculate and save fjord-integrated values
s.IT = p.W*p.L*trapz(f.zi,s.I); % fjord iceberg volume
s.MT = p.W*p.L*trapz(f.zi,s.M.*s.I); % total iceberg melt flux
s.ET = p.W*p.L*trapz(f.zi,p.E0*s.I); % total iceberg export flux

% return forcing on same timestepping
f.Ss = f.Ss(:,1:int:end-1);
f.Ts = f.Ts(:,1:int:end-1);
f.Qsg = f.Qsg(1:int:end-1);
f.D = f.D(1:int:end-1);

%% Save output
% save(['./output_',name,'/out_boxmodel.mat'],'s','f');
