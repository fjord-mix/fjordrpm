function [s,f] = boxmodel_v4(p,f,a,t);

% INPUTS
% p - constant parameters structure
% f - forcings structure
% a - initial conditions structure
% t - time variable

% OUTPUTS
% s - solution stucture, contents below

% initialise variables
H(:,1) = a.H0; % thickness
T(:,1) = a.T0; % temperature
S(:,1) = a.S0; % salinity
I(:,1) = a.I0; % iceberg concentration
V(:,1) = H(:,1).*double(p.W).*double(p.L); % volume
VT(:,1) = V(:,1).*T(:,1); % heat
VS(:,1) = V(:,1).*S(:,1); % salt

if p.plot_runtime
    hf_track = monitor_boxmodel([],1,H,T,S,f);
end
% the main loop
for i=1:length(t)-1,    

    % calculate plume fluxes
    [QVg(:,i),QTg(:,i),QSg(:,i)] = ...
        plume_fluxes(H(:,i),T(:,i),S(:,i),f.Qsg(i),p);
    if ~isreal(H(:,i))
    disp(['it became imaginary at timestep ', num2str(i)])
    end
    % calculate shelf fluxes
    [QVs(:,i),QTs(:,i),QSs(:,i),Se(:,i),Te(:,i),phi(:,i)] = ...
        shelf_fluxes(H(:,i),T(:,i),S(:,i),f.zs,f.Ts(:,i),f.Ss(:,i),f.Qsg(i),p);

    % calculate vertical mixing fluxes
    [QVk(:,i),QTk(:,i),QSk(:,i)] = ...
        mixing_fluxes(H(:,i),T(:,i),S(:,i),QVg(:,i),QVs(:,i),p,i);
    
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
        hf_track = monitor_boxmodel(hf_track,i,H,T,S,f);
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

%% function to calculate plume fluxes
function [QpV0,QpT0,QpS0] = plume_fluxes(H0,T0,S0,Qsg0,p);

    if Qsg0==0 | p.P0==0, % i.e. if no plume
    
        QpV0 = 0*H0;
        QpT0 = 0*H0;
        QpS0 = 0*H0;    
    
    else
    
        % find box model layer containing grounding line
        ints = cumsum(H0);
        kgl = min(find(ints>=abs(p.zgl)-1e-6));
        if isempty(kgl), kgl = 4; end % in case the gl is at the bottom of the fjord
    
        % plume properties at grounding line
        k = kgl;
        Qp(k) = Qsg0;
        Sp(k) = 0;
        Tp(k) = 0;
        gp(k) = p.g*(p.betaS*(S0(k)-Sp(k))-p.betaT*(T0(k)-Tp(k)));
    
        % properties at first interface about grounding line
        % need special treatment because box might be partial if
        % grounding line does not coincide with box boundaries
        k = k-1;
        % if k < 1, k=1; end % provisory fix to avoid the "interface" being at a box that does not exist        
        try
        Qp(k) = Qp(k+1) + p.P0^(2/3)*Qp(k+1)^(1/3)*gp(k+1)^(1/3)*(abs(p.zgl)-ints(k));
        Tp(k) = (Qp(k+1)*Tp(k+1)+(Qp(k)-Qp(k+1))*T0(k+1))/Qp(k);
        Sp(k) = (Qp(k+1)*Sp(k+1)+(Qp(k)-Qp(k+1))*S0(k+1))/Qp(k);
        gp(k) = p.g*(p.betaS*(S0(k)-Sp(k))-p.betaT*(T0(k)-Tp(k)));
        catch
            disp('there is a problem here')
        end
        % apply to successive interfaces higher provided plume is still rising
        while gp(k)>0 & k>1,
            k = k-1;
            Qp(k) = Qp(k+1) + p.P0^(2/3)*Qp(k+1)^(1/3)*gp(k+1)^(1/3)*H0(k+1);
            Tp(k) = (Qp(k+1)*Tp(k+1)+(Qp(k)-Qp(k+1))*T0(k+1))/Qp(k);
            Sp(k) = (Qp(k+1)*Sp(k+1)+(Qp(k)-Qp(k+1))*S0(k+1))/Qp(k);
            gp(k) = p.g*(p.betaS*(S0(k)-Sp(k))-p.betaT*(T0(k)-Tp(k)));
        end
    
        % now calculate the resulting box fluxes
        knb = find(gp<0); if isempty(knb), knb=1; end % find neutral buoyancy box
        % flux in boxes below grounding line and above neutral buoyancy are 0
        inds = find([1:length(H0)]>kgl | [1:length(H0)]<knb);        
        QpV0(inds) = 0;
        QpT0(inds) = 0;
        QpS0(inds) = 0;
        
        % flux in boxes from grounding line to below neutral buoyancy
        for k=kgl:-1:knb+1,
            QpV0(k) = Qp(k)-Qp(k-1);
            QpT0(k) = QpV0(k)*T0(k);
            QpS0(k) = QpV0(k)*S0(k);
        end
        % flux into neutral buoyancy box
        QpV0(knb) = Qp(knb);
        QpT0(knb) = Qp(knb)*Tp(knb);
        QpS0(knb) = Qp(knb)*Sp(knb);
    
    end

end

%% function to calculate shelf fluxes
function [QVs0,QTs0,QSs0,Se0,Te0,phi0] = shelf_fluxes(H0,T0,S0,zs,Ts,Ss,Qsg0,p),

    if p.C0==0,
        
        QVs0 = 0*H0;
        QTs0 = 0*H0;
        QSs0 = 0*H0;
        Se0 = NaN*H0;
        Te0 = NaN*H0;
        phi0 = NaN*H0;
    
    else
    
        % calculate mean shelf T/S over box model layers
        ints = [0;-cumsum(H0)];
        zs0 = unique(sort([zs,-cumsum(H0)']));        
        Ss0 = interp1(zs,Ss,zs0,'pchip','extrap');
        Ts0 = interp1(zs,Ts,zs0,'pchip','extrap');
        
        for k=1:length(ints)-1,
            inds = find(zs0<=ints(k) & zs0>=ints(k+1));
            if length(inds) == 1 % if there is only one data point, no need to average it
                Se0(k) = Ss0(inds);
                Te0(k) = Ts0(inds);
            else
                Se0(k) = trapz(zs0(inds),Ss0(inds))/H0(k);
                Te0(k) = trapz(zs0(inds),Ts0(inds))/H0(k);
            end
        end
        Se0 = Se0';
        Te0 = Te0';
    
        % get fjord to shelf reduced gravity
        for k=1:length(ints)-1,
            gp(k) = p.g*(p.betaS*(S0(k)-Se0(k))-p.betaT*(T0(k)-Te0(k)));
        end
    
        % calculate potentials
        for k=1:3,
            if k==1,
                phi0(k) = gp(k)*H0(k)/2;
            else
                phi0(k) = phi0(k-1)+gp(k-1)*H0(k-1)/2+gp(k)*H0(k)/2;
            end
        end
    
        % fluxes before barotropic compensation
        Q = p.C0*p.W*H0(1:3).*phi0'/p.L;
    
        % fluxes after ensuring depth mean = QSg0
        QVs0 = Q + H0(1:3)*(Qsg0-sum(Q))/sum(H0(1:3));
        QVs0(4) = 0;
    
        % resulting heat/salt fluxes
        QTs0 = (QVs0>0).*QVs0.*T0 + (QVs0<0).*QVs0.*Te0;
        QSs0 = (QVs0>0).*QVs0.*S0 + (QVs0<0).*QVs0.*Se0;
    
    end

end

%% function to calculate vertical mixing fluxes
function [QVk0,QTk0,QSk0] = mixing_fluxes(H0,T0,S0,QVg0,QVs0,p,i_model)
if i_model==22
    disp('w will become imaginary here')
end

    if p.K0==0,
    
        QVk0 = 0*H0;
        QTk0 = 0*H0;
        QSk0 = 0*H0;
    
    else
       
        % vector for doing numerical integral
        x = linspace(0,p.L,10);
    
        % loop over 1-2 interface and 2-3 interface
        for k=1:2,
            % velocity in lower box
            ul = (QVg0(k+1)+(QVs0(k+1)-QVg0(k+1))*x/p.L)/(p.W*H0(k+1));
            % velocity in higher box
            uu = (QVg0(k)+(QVs0(k)-QVg0(k))*x/p.L)/(p.W*H0(k));
            % buoyancy jump between boxes
            B = p.g*(p.betaS*(S0(k+1)-S0(k))-p.betaT*(T0(k+1)-T0(k)));
            % richardson number
            R = 0.5*(H0(k+1)+H0(k))*B./(ul-uu).^2;
            % vertical entrainment velocity
            % capped at 4e-5
            w = min(p.K0*abs(ul-uu).*R.^(-0.75),p.wmax);
            % fluxes
            Q(k) = sign(mean(abs(uu)-abs(ul)))*p.W*trapz(x,w);
            QT(k) = (Q(k)>0)*Q(k)*T0(k+1) + (Q(k)<0)*Q(k)*T0(k);
            QS(k) = (Q(k)>0)*Q(k)*S0(k+1) + (Q(k)<0)*Q(k)*S0(k);
        end
    
        % final fluxes
        QVk0 = [Q(1);-Q(1)+Q(2);-Q(2);0];
        QTk0 = [QT(1);-QT(1)+QT(2);-QT(2);0];
        QSk0 = [QS(1);-QS(1)+QS(2);-QS(2);0];
    
    end
       
end

%% function to calculate "artificial" fluxes
function [QVb0,QTb0,QSb0] = artificial_fluxes(QT0,H0,V0,T0,S0,zs,Ts,Ss,p);

    % maintain volume of box 4
    QV43 = QT0(4);
    QT43 = (QV43>0)*QV43*T0(4) + (QV43<0)*QV43*T0(3);
    QS43 = (QV43>0)*QV43*S0(4) + (QV43<0)*QV43*S0(3);

    % resulting box-to-box vectors
    QVb0 = [0;0;QV43;-QV43];
    QTb0 = [0;0;QT43;-QT43];
    QSb0 = [0;0;QS43;-QS43];

    % do layer nudging followed by minimum thickness
    if ~isnan(p.trelax), % if layer nudging active

        % get depth of layer boundaries
        Z0 = cumsum(H0);

        % get desired depth based on salinity boundaries on shelf
        [~,minind] = min(abs(Ss-p.S12)); z12 = abs(zs(minind));
        [~,minind] = min(abs(Ss-p.S23)); z23 = abs(zs(minind));

        % nudging fluxes
        QV21 = p.W*p.L*(z12-Z0(1))/p.trelax;
        QV32 = p.W*p.L*(z23-Z0(2))/p.trelax;

        % associated heat and salt fluxes
        QT21 = (QV21>0)*QV21*T0(2) + (QV21<0)*QV21*T0(1);
        QS21 = (QV21>0)*QV21*S0(2) + (QV21<0)*QV21*S0(1);
        QT32 = (QV32>0)*QV32*T0(3) + (QV32<0)*QV32*T0(2);
        QS32 = (QV32>0)*QV32*S0(3) + (QV32<0)*QV32*S0(2);

        % update box-to-box vectors
        QVb0 = QVb0 + [QV21;QV32-QV21;-QV32;0];
        QTb0 = QTb0 + [QT21;QT32-QT21;-QT32;0];
        QSb0 = QSb0 + [QS21;QS32-QS21;-QS32;0];        

    end

    if ~isnan(p.Hmin), % if minimum thickness active

        % update the total fluxes vector for possible nudging
        QT0 = QT0 + QVb0;
    
        % reset 21 and 32 variables
        QV21 = 0;
        QT21 = 0;
        QS21 = 0;
        QV32 = 0;
        QT32 = 0;
        QS32 = 0;
    
        % routine to avoid going below minimum thickness
        % this is awkward but the following should cover all possibilities
        Vtend = p.dt*p.sid*QT0;
        inds = find(V0+Vtend<p.W*p.L*p.Hmin);
        
        % case where only box 2 is too thin
        if isequal(inds,2),
            % if box 1 is bigger than 3, take from 1
            if V0(1)>V0(3),
                QV21 = QV21-(p.W*p.L*p.Hmin-V0(2)-Vtend(2))/(p.dt*p.sid);
                QT21 = QV21*T0(1);
                QS21 = QV21*S0(1);
            % else take from 3
            else    
                QV32 = QV32+(p.W*p.L*p.Hmin-V0(2)-Vtend(2))/(p.dt*p.sid);
                QT32 = QV32*T0(3);
                QS32 = QV32*S0(3);
            end
        
        % case where 1 and 3 are too thin
        % take from box 2
        elseif isequal(inds,[1,3]),
            % box 2 to 1
            QV21 = QV21+(p.W*p.L*p.Hmin-V0(1)-Vtend(1))/(p.dt*p.sid);
            QT21 = QV21*T0(2);
            QS21 = QV21*S0(2);
            % box 2 to 3
            QV32 = QV32-(p.W*p.L*p.Hmin-V0(3)-Vtend(3))/(p.dt*p.sid);
            QT32 = QV32*T0(2);
            QS32 = QV32*S0(2);
        
        % case where 1 is too thin (and possibly 2, but not 3)
        elseif any(inds==1),
            % box 2 to 1
            QV21 = QV21+(p.W*p.L*p.Hmin-V0(1)-Vtend(1))/(p.dt*p.sid);
            QT21 = QV21*T0(2);
            QS21 = QV21*S0(2);
            % recalculate tendency for box 2
            Vtend(2) = p.dt*p.sid*(QT0(2)+QVb0(2)-QV21);
            % if tendency would make box 2 too thin, take from 3
            if V0(2)+Vtend(2)<p.W*p.L*p.Hmin,
                QV32 = QV32+(p.W*p.L*p.Hmin-V0(2)-Vtend(2))/(p.dt*p.sid);
                QT32 = QV32*T0(3);
                QS32 = QV32*S0(3);
            end
    
        % case where 3 is too thin (and possibly 2, but not 1)
        elseif any(inds==3),
            % box 2 to 3
            QV32 = QV32-(p.W*p.L*p.Hmin-V0(3)-Vtend(3))/(p.dt*p.sid);
            QT32 = QV32*T0(2);
            QS32 = QV32*S0(2);        
             % recalculate tendency for box 2
            Vtend(2) = p.dt*p.sid*(QT0(2)+QVb0(2)+QV32);   
            % if tendency would make box 2 too thin, take from 1
            if V0(2)+Vtend(2)<p.W*p.L*p.Hmin,
                QV21 = QV21-(p.W*p.L*p.Hmin-V0(2)-Vtend(2))/(p.dt*p.sid);
                QT21 = QV21*T0(1);
                QS21 = QV21*S0(1);
            end
    
        end
    
        % update box-to-box vector
        QVb0 = QVb0 + [QV21;QV32-QV21;-QV32;0];
        QTb0 = QTb0 + [QT21;QT32-QT21;-QT32;0];
        QSb0 = QSb0 + [QS21;QS32-QS21;-QS32;0];

    end

end

%% function to calculate iceberg fluxes
function [QIi0,QTi0,QSi0,M0] = iceberg_fluxes(H0,T0,S0,I0,zi,p);

    if p.M0==0,

        QIi0 = 0*H0;
        QTi0 = 0*H0;
        QSi0 = 0*H0;
        M0 = 0*zi;

    else

        % melting in boxes
        ints = [0;cumsum(H0)];
        zj = 0.5*(ints(1:end-1)+ints(2:end)); % mean depth of box
        Tf = p.l1*S0 + p.l2 + p.l3*zj; % local freezing point
        % get vector of iceberg concentration that resolves box boundaries
        zi0 = unique(sort([zi,-cumsum(H0)']));
        Ii0 = interp1(zi,I0,zi0,'pchip','extrap');
        % do numerical integral to get melt fluxes
        for k=1:length(ints)-1,
            inds = find(zi0<=ints(k) & zi0>=ints(k+1));
            QIi0(k) = p.W*p.L*p.M0*(T0(k)-Tf(k))*trapz(zi0(inds),Ii0(inds)); % modified Tf0 -> Tf
        end 
        QTi0 = -QIi0*p.l/p.cw;
        QSi0 = -QIi0.*S0'; % modified S0 -> S0'

        % melt profile defined on zi reference depths
        M0 = 0*zi;
        for k=1:length(ints)-1,
            inds = find(zi<=ints(k) & zi>=ints(k+1));
            M0(inds) = p.M0*(T0(k)-Tf(k));
        end

    end

end

