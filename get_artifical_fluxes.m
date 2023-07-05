%% function to calculate "artificial" fluxes
function [QVb0,QTb0,QSb0] = artificial_fluxes(QT0,H0,V0,T0,S0,zs,Ts,Ss,p)

    % maintain volume of below-sill box
    if p.sill==1,

        QVsill = QT0(p.N+p.sill);
        QTsill = (QVsill>0)*QVsill*T0(p.N+p.sill) + (QVsill<0)*QVsill*T0(p.N);
        QSsill = (QVsill>0)*QVsill*S0(p.N+p.sill) + (QVsill<0)*QVsill*S0(p.N);
    
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

        % get depth of layer boundaries
        Z0 = cumsum(H0);

        % get desired depth based on salinity boundaries on shelf
        for k=1:p.N-1,

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
        if p.sill==1,
            QVnudge = [QVnudge,0];
            QTnudge = [QTnudge,0];
            QSnudge = [QSnudge,0];
        end
        QVb0 = QVb0 + [QVnudge,0]'-[0,QVnudge]';    
        QTb0 = QTb0 + [QTnudge,0]'-[0,QTnudge]'; 
        QSb0 = QSb0 + [QSnudge,0]'-[0,QSnudge]'; 

    end

%     if ~isnan(p.Hmin) % if minimum thickness active
% 
%         % update the total fluxes vector for possible nudging
%         QT0 = QT0 + QVb0;
%     
%         % routine to avoid going below minimum thickness
%         Vtend = p.dt*p.sid*QT0;
%         mininds = find(V0+Vtend<p.W*p.L*p.Hmin);
%         inds = find(V0+Vtend>=p.W*p.L*p.Hmin);
%         QVmin = zeros(1,p.N+p.sill-1);
% 
%         % for each box about to go below min thickness, find the nearest
%         % box that is not about to go below min thickness and take volume
%         % from that, which may require going through an intermediate box
%         % this is awkward to code so leave for now and put
%         % catch into the code to stop if box thickness goes below 0
% 
%     end

end