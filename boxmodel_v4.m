function [s,f] = boxmodel_v4(p,f,a,t)

% INPUTS
% p - constant parameters structure
% f - forcings structure
% a - initial conditions structure
% t - time variable

% OUTPUTS
% s - solution stucture, contents below
s.status=0; % initial assumption that the model ran successfully

% initialise variables
H(:,1) = a.H0; % thickness
T(:,1) = a.T0; % temperature
S(:,1) = a.S0; % salinity
I(:,1) = a.I0; % iceberg concentration
V(:,1) = H(:,1).*double(p.W).*double(p.L); % volume
VT(:,1) = V(:,1).*T(:,1); % heat
VS(:,1) = V(:,1).*S(:,1); % salt

if p.plot_runtime
    % hf_track = monitor_boxmodel([],1,H,T,S,f);
    hf_track = show_box_model([],1,t,H,T,S,[],[],[],[],f);
end
% the main loop
for i=1:length(t)-1    
    try

    % calculate plume fluxes
    [QVg(:,i),QTg(:,i),QSg(:,i)] = ...
        plume_fluxes(H(:,i),T(:,i),S(:,i),f.Qsg(i),p,i);
    
    % calculate shelf fluxes
    [QVs(:,i),QTs(:,i),QSs(:,i),Se(:,i),Te(:,i),phi(:,i)] = ...
        shelf_fluxes(H(:,i),T(:,i),S(:,i),f.zs,f.Ts(:,i),f.Ss(:,i),f.Qsg(i),p);

    % calculate vertical mixing fluxes
    [QVk(:,i),QTk(:,i),QSk(:,i)] = ...
        mixing_fluxes(H(:,i),T(:,i),S(:,i),QVg(:,i),QVs(:,i),p);
    
    % calculate "artificial" fluxes
    [QVb(:,i),QTb(:,i),QSb(:,i)] = ...
        artificial_fluxes(QVg(:,i)-QVs(:,i)+QVk(:,i),H(:,i),V(:,i),T(:,i),S(:,i),f.zs,f.Ts(:,i),f.Ss(:,i),p);

    % calculate iceberg fluxes
    [QIi(:,i),QTi(:,i),QSi(:,i),M(:,i)] = ...
        iceberg_fluxes(H(:,i),T(:,i),S(:,i),I(:,i),f.zi,p);
    

    % step fjord forwards
    dt = t(i+1)-t(i);
    V(:,i+1)  = V(:,i)  + dt*p.sid*(QVg(:,i)-QVs(:,i)+QVk(:,i)+QVb(:,i));
    VT(:,i+1) = VT(:,i) + dt*p.sid*(QTg(:,i)-QTs(:,i)+QTk(:,i)+QTb(:,i)+QTi(:,i));
    VS(:,i+1) = VS(:,i) + dt*p.sid*(QSg(:,i)-QSs(:,i)+QSk(:,i)+QSb(:,i)+QSi(:,i));

    % calculate thicknesses and tracers
    H(:,i+1) = V(:,i+1)/(p.W*p.L);
    T(:,i+1) = VT(:,i+1)./V(:,i+1);
    S(:,i+1) = VS(:,i+1)./V(:,i+1);
    
    % step icebergs forwards
    I(:,i+1) = I(:,i) + dt*p.sid*((f.D(i)/(p.W*p.L))*f.xi-M(:,i).*I(:,i)-p.E0*I(:,i));

    % plot model evolution (mainly debugging)
    if p.plot_runtime
        % hf_track = monitor_boxmodel(hf_track,i,H,T,S,f);
        hf_track = show_box_model(hf_track,i,t,H,T,S,QVs,QVg,QVk,QVb,f);
    end

        % Check for possible blow up
    catch
        disp(['Model error at timestep ', num2str(i)])
        disp('Breaking out of loop and saving outputs')
        s.status=1; % status == 1 means there was an error
        break
    end
end

% put solution into output structure
% but just save ~daily values
% as otherwise high time resolution results in large output files
dtdaily = 1;
int = round(dtdaily/dt);
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
% derive 21 and 32 fluxes from these
s.QV21 = s.QVk(1,:);
s.QV32 = s.QVk(1,:)+s.QVk(2,:);

% artificial fluxes
s.QVb = QVb(:,1:int:end);
s.QTb = QTb(:,1:int:end);
s.QSb = QSb(:,1:int:end);
% 43 flux
s.QV43 = -s.QVb(4,:);

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

end

