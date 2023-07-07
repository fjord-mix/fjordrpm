function [QVb0,QTb0,QSb0] = get_artificial_fluxes(QV0,H0,V0,T0,S0,zs,Ts,Ss,p)

% GET_ARTIFICIAL_FLUXES Compute artificial fluxes.
%   [QVB0,QTB0,QSB0] = GET_ARTIFICIAL_FLUXES(QV0,H0,V0,T0,S0,ZS,TS,SS,P)
%   computes the artifcial fluxes for the given parameters. 

if p.sill == 1
    QVsill = QV0(p.N+p.sill);
    QTsill = (QVsill>0)*QVsill*T0(p.N+p.sill)+(QVsill<0)*QVsill*T0(p.N);
    QSsill = (QVsill>0)*QVsill*S0(p.N+p.sill)+(QVsill<0)*QVsill*S0(p.N);
    % resulting box-to-box vectors
    QVb0 = [zeros(p.N-p.sill,1);QVsill;-QVsill];
    QTb0 = [zeros(p.N-p.sill,1);QTsill;-QTsill];
    QSb0 = [zeros(p.N-p.sill,1);QSsill;-QSsill];
else
    QVb0 = 0*H0;
    QTb0 = 0*H0;
    QSb0 = 0*H0;
end

% do layer nudging followed by minimum thickness
if ~isnan(p.trelax) % if layer nudging active
    % Preallocate variables
    [zn, QVnudge, QTnudge, QSnudge] = deal(zeros(1, p.N-1));
    % get depth of layer boundaries
    Z0 = cumsum(H0);
    % get desired depth based on salinity boundaries on shelf
    for k = 1:p.N-1
        % desired depth
        [~,minind] = min(abs(Ss-p.Snudge(k)));
        zn(k) = abs(zs(minind));
        % nudging flux
        QVnudge(k) = p.W*p.L*(zn(k)-Z0(k))/(p.trelax*p.sid);
        % associated heat and salt fluxes
        QTnudge(k) = (QVnudge(k)>0)*QVnudge(k)*T0(k+1)+(QVnudge(k)<0)*QVnudge(k)*T0(k);
        QSnudge(k) = (QVnudge(k)>0)*QVnudge(k)*S0(k+1)+(QVnudge(k)<0)*QVnudge(k)*S0(k);
    end

    % update box-to-box vectors
    if p.sill == 1
        QVnudge = [QVnudge,0];
        QTnudge = [QTnudge,0];
        QSnudge = [QSnudge,0];
    end
    QVb0 = QVb0+[QVnudge,0]'-[0,QVnudge]';
    QTb0 = QTb0+[QTnudge,0]'-[0,QTnudge]';
    QSb0 = QSb0+[QSnudge,0]'-[0,QSnudge]';

end

if ~isnan(p.Hmin) % if minimum thickness active

    % update the total fluxes vector for possible nudging
    QV0 = QV0 + QVb0;

    % routine to avoid going below minimum thickness
    Vtend = p.dt*p.sid*QV0;
    mininds = find(V0+Vtend<p.W*p.L*p.Hmin);    
    QVmin = zeros(p.N+p.sill,1);            

    % for each box about to go below min thickness, find the nearest
    % (and thickest) box that is not about to go below min thickness 
    % and take volume from that        
    while ~isempty(mininds)
        Has = H0(1:p.N); % thicknesses of layers above the sill
        i_from_box = zeros(size(QVmin)); % we will need to keep track of which box we take the water from for nudging
        for i=1:length(mininds)
            i_thin = mininds(i);

            % will take water from the thickest adjacent layer, unless
            % there is none on top, or if the one below is the sill
            if i_thin == p.N
                i_thickest = i_thin-1;
            elseif i_thin ==1
                i_thickest = i_thin+1;
            else
                thick_layer = max(Has(i_thin-1),Has(i_thin+1));
                i_thickest = find(Has == thick_layer);
            end
            i_from_box(i_thin) = i_thickest; % recording from which box the water was taken
            
            QVmin(i_thin)     = QVmin(i_thin)    +(p.W*p.L*p.Hmin-V0(i_thin)-Vtend(i_thin))/(p.dt*p.sid);
            QVmin(i_thickest) = QVmin(i_thickest)-(p.W*p.L*p.Hmin-V0(i_thin)-Vtend(i_thin))/(p.dt*p.sid);                
        end
        Vtend = p.dt*p.sid*(QV0+QVb0+QVmin);
        mininds = find(V0+Vtend<p.W*p.L*p.Hmin);
    end

    if exist('i_from_box','Var') && ~isempty(find(QVmin ~=0,1))
        QVb0 = QVb0 + QVmin;
        for k=1:p.N
            i_from = i_from_box(k);
            if i_from ~= 0
                QTb0(k) = QTb0(k) + (QVmin(k)>0)*QVmin(k)*T0(i_from)+(QVmin(k)<0)*QVmin(k)*T0(k);
                QSb0(k) = QSb0(k) + (QVmin(k)>0)*QVmin(k)*S0(i_from)+(QVmin(k)<0)*QVmin(k)*S0(k);        
            end
        end 
    end

end

end