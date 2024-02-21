function [s,f] = boxmodel(p,t,f,a,path_out)

% BOXMODEL Box model simulation.
%   [S,F] = BOXMODEL(P,T,F,A,PATH_OUT) runs the box model simulation for parameters structure P,
%   time T, forcings structure F, initial conditions structure A. Returns solution structure
%   S and forcing structure F in the same time steps as S. If PATH_OUT is
%   specified, will save a file.
%   If F and/or A are not specified, default idealised values will be used.

s.status = 0; % initial assumption that the model ran successfully

%% Initialise variables
if nargin < 3, f = get_idealised_forcing(p, t); end
if isempty(f), f = get_idealised_forcing(p, t); end % we cannot use an OR statement here

if nargin < 4, a = get_initial_conditions(p, f); end
if isempty(a), a = get_initial_conditions(p, f); end % we cannot use an OR statement here
dt = p.dt;

H(:,1) = a.H0; % thickness
T(:,1) = a.T0; % temperature
S(:,1) = a.S0; % salinity
I(:,1) = a.I0; % iceberg concentration
V(:,1) = H(:,1).*double(p.W).*double(p.L); % volume
VT(:,1) = V(:,1).*T(:,1); % heat
VS(:,1) = V(:,1).*S(:,1); % salt

% Preallocate variables
[QVg,QTg,QSg,QVs,QTs,QSs,QVk,QTk,QSk,QVmi,QTmi,QSmi,QVb,QTb,QSb,QIi,QTi,QSi,Te,Se] = deal(zeros(size(H,1),length(t)-1));
phi = zeros(size(H,1)-p.sill,length(t)-1);
M = zeros(length(f.zi),length(t)-1);
homogenisation_timestep = zeros(1, length(t)-1);

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
if p.sill == 1 && (H(end,1) - (p.H-abs(p.silldepth))) > 1e-4
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

%%
if p.plot_runtime
    % hf_track = monitor_boxmodel([],1,H,T,S,f);
    % hf_track = show_boxmodel([],1,t,H,T,S,[],[],[],[],f);
    s_bnds = [min(f.Ss(:)) max(f.Ss(:))+0.1];
    plot_debug_profile(0,t,f,p,H,S,s_bnds);
end

%% The main loop
for i = 1:length(t)-1
    
    % check if the fjord stratification is unstable, and homogenise layer
    % properties if it is
    for k=1:p.N-1        
        % buoyancy jump between boxes
        B = p.g*(p.betaS*(S(k+1,i)-S(k,i))-p.betaT*(T(k+1,i)-T(k,i)));        
        if B < 0
            [T(:,i),S(:,i),V(:,i),VT(:,i),VS(:,i),H(:,i)] ...
            = homogenise_layers(V(:,i),T(:,i),S(:,i),VT(:,i),VS(:,i),[k,k+1],p.L,p.W);
        end
    end
    % if there is a sill
    if p.sill == 1
        % check the buoyancy jump between the below sill layer and layer N
        k = p.N;
        B = p.g*(p.betaS*(S(k+1,i)-S(k,i))-p.betaT*(T(k+1,i)-T(k,i)));
        % if the buoyancy jump is negative, homogenise the heat and salt
        % properties of the two layers but don't change their volume 
        if B < 0
            inds = [k, k+1];
            T(inds,i) = sum(T(inds,i).*V(inds,i))./sum(V(inds,i));
            S(inds,i) = sum(S(inds,i).*V(inds,i))./sum(V(inds,i));
            % Recompute heat and salt content 
            VT(inds,i) = V(inds,i).*T(inds,i);
            VS(inds,i) = V(inds,i).*S(inds,i);
        end
    end
    
    counter = 0;
    while counter <= p.N
        % Compute fluxes. Check to see if any layer has collapsed, and if so homogenise and
        % recompute fluxes.

        % Compute fluxes
        [...
            QVg(:,i),QTg(:,i),QSg(:,i),...
            QVs(:,i),QTs(:,i),QSs(:,i),Se(:,i),Te(:,i),phi(:,i),...
            QVk(:,i),QTk(:,i),QSk(:,i),...
            QVb(:,i),QTb(:,i),QSb(:,i),...
            QIi(:,i),QTi(:,i),QSi(:,i),M(:,i),QVmi(:,i),QTmi(:,i),QSmi(:,i)] ...
            = compute_fluxes(...
            H(:,i),T(:,i),S(:,i),f.Qsg(i),p,f.zs,f.Ts(:,i),f.Ss(:,i), ...
            V(:,i),I(:,i),f.zi);
        
        [homogenisation_flag, H(:,i),T(:,i),S(:,i),V(:,i), VT(:,i),VS(:,i)] = ...
            homogenise_thin_layers(V(:,i),T(:,i),S(:,i), VT(:,i), VS(:,i), H(:,i), dt, p.sid, QVg(:,i), ...
            QVs(:,i), QVk(:,i), QVb(:,i), p.W, p.L, p.N, p);

        if homogenisation_flag == false
            % If no homogenisation has occurred, timestep forwards
            break
        else
            % If homogenisation has occurred, recompute fluxes and go through
            % the loop again.
            counter = counter +1;
        end
    end
    if counter > p.N
        % If p.N iterations of homogenisation have occurred and layers are
        % still collapsing, break because the model is unstable
        error("Error: collapsing layers still exist after homogenisation, model is unstable")
    end

    % Store timestep at which this has occured and how many iterations were
    % needed
    s.homogenisation_timestep(i) = counter;

    % Step fjord forwards.
    
    % dt = t(i+1)-t(i); % replaced by pre-defined dt because of problems when running in parallel
    V(:,i+1)  = V(:,i)+dt*p.sid*(QVg(:,i)-QVs(:,i)+QVk(:,i)+QVmi(:,i)+QVb(:,i));
    VT(:,i+1) = VT(:,i)+dt*p.sid*(QTg(:,i)-QTs(:,i)+QTk(:,i)+QTi(:,i)+QTmi(:,i)+QTb(:,i));
    VS(:,i+1) = VS(:,i)+dt*p.sid*(QSg(:,i)-QSs(:,i)+QSk(:,i)+QSi(:,i)+QSmi(:,i)+QSb(:,i));

    % Calculate thicknesses and tracers.
    H(:,i+1) = V(:,i+1)/(p.W*p.L);
    T(:,i+1) = VT(:,i+1)./V(:,i+1);
    S(:,i+1) = VS(:,i+1)./V(:,i+1);
    % if ~isempty(find(S(:,end) > 35,1))
    %     disp('started becoming unstable')
    % end

    % Step icebergs forwards.
    I(:,i+1) = I(:,i);%+dt*p.sid*((f.D(i)/(p.W*p.L))*f.xi-M(:,i).*I(:,i)-p.E0*I(:,i));

    % real-time nudging, i.e., nudging values are updated to mimic the current shelf conditions
    if ~isnan(p.trelax) && p.real_time_nudge
        p.Snudge = get_interface_salinities(f.zs,f.Ts(:,i),f.Ss(:,i),p);
    end

    % Plot model evolution (mainly debugging).
    if p.plot_runtime
        % hf_track = monitor_boxmodel(hf_track,i,H,T,S,f);
        % hf_track = show_boxmodel([],i,t,H,T,S,QVs,QVg,QVk,QVb,f);
        plot_debug_profile(i,t,f,p,H,S,[]);
    end

    % if ~isempty(find(H(:,i+1) < p.Hmin,1))
    %     thinnest_layer=find(H(:,i+1) < p.Hmin,1);
    %     fprintf('Error: layer %d thickness dropped below p.Hmin at time step %d\n',thinnest_layer,i);
    %     s.status = 1; % status == 1 means there was an error
    %     if thinnest_layer > 2
    %         disp('we want to save this')
    %     end
    %     break
    % end

    % Break from loop if the sill layer is not acting as it was supposed to
    % Check bottom box is consistent with sill depth.
    if p.sill == 1 && (H(end,i+1) - (p.H-abs(p.silldepth))) > 1e-4
        disp('Error: when p.sill=1, bottom box must have thickness p.H-p.silldepth');
        s.status = 1; % status == 1 means there was an error
        break
    end
    % Check sum of layer thicknesses is equal to fjord depth.
    if abs(sum(H(:,i+1))-p.H) > 1e-10
        disp('Error: box thicknesses must sum to fjord depth');
        s.status = 1; % status == 1 means there was an error
        break
    end
    % check if something starts to become unstable
    % if ~isempty(find(S(:,i+1) > 45,1))
    %     [sind,~]=find(S(:,i+1) > 38);
    %     fprintf('Salinity went above 38 in layer %d at time step %d\n',sind,i)
    % end
    % if ~isempty(find(S(:,i+1) <0,1))
    %     [sind,~]=find(S(:,i+1) <0,1);
    %     fprintf('Salinity went negative in layer %d\n',sind)
    % end

end

%% Put solution into output structure
% but just save ~daily values
% as otherwise high time resolution results in large output files
dtdaily = max(1, p.dt);
if s.status == 0
    int = round(dtdaily/dt);
else
    int = 1; % if something went wrong, we want all time steps to properly understand what happened
end
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

% artificial fluxes
s.QVb = QVb(:,1:int:end);
s.QTb = QTb(:,1:int:end);
s.QSb = QSb(:,1:int:end);

% iceberg fluxes
s.QIi = QIi(:,1:int:end);
s.QTi = QTi(:,1:int:end);
s.QSi = QSi(:,1:int:end);
s.QVmi = QVmi(:, 1:int:end);
s.M = M(:,1:int:end);

% for iceberg fluxes also calculate and save fjord-integrated values
s.IT = p.W*p.L*trapz(f.zi,s.I); % fjord iceberg volume
s.MT = p.W*p.L*trapz(f.zi,s.M(:,1:size(s.I,2)).*s.I); % total iceberg melt flux
s.ET = p.W*p.L*trapz(f.zi,p.E0*s.I); % total iceberg export flux

% return forcing on same time step as forcings (in results structure to prevent overwriting)
s.Ss = f.Ss(:,1:int:end-1);
s.Ts = f.Ts(:,1:int:end-1);
s.Qsg = f.Qsg(1:int:end-1);
s.D = f.D(1:int:end-1);

%% Save output if a path+file name are provided
if nargin > 4
    fjord_output.s = s;
    fjord_output.f = f;
    fjord_output.t = t;
    fjord_output.p = p;
    fjord_output.a = a;
    save(path_out,'fjord_output','-v7.3'); % v7.3 allows large files (> 2GB), which might happen in very long runs
end

end
