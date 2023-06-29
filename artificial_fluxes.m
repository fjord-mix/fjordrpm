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

    if ~isnan(p.Hmin) % if minimum thickness active

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
        if isequal(inds,2)
            % if box 1 is bigger than 3, take from 1
            if V0(1)>V0(3)
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
        elseif isequal(inds,[1,3])
            % box 2 to 1
            QV21 = QV21+(p.W*p.L*p.Hmin-V0(1)-Vtend(1))/(p.dt*p.sid);
            QT21 = QV21*T0(2);
            QS21 = QV21*S0(2);
            % box 2 to 3
            QV32 = QV32-(p.W*p.L*p.Hmin-V0(3)-Vtend(3))/(p.dt*p.sid);
            QT32 = QV32*T0(2);
            QS32 = QV32*S0(2);
        
        % case where 1 is too thin (and possibly 2, but not 3)
        elseif any(inds==1)
            % box 2 to 1
            QV21 = QV21+(p.W*p.L*p.Hmin-V0(1)-Vtend(1))/(p.dt*p.sid);
            QT21 = QV21*T0(2);
            QS21 = QV21*S0(2);
            % recalculate tendency for box 2
            Vtend(2) = p.dt*p.sid*(QT0(2)+QVb0(2)-QV21);
            % if tendency would make box 2 too thin, take from 3
            if V0(2)+Vtend(2)<p.W*p.L*p.Hmin
                QV32 = QV32+(p.W*p.L*p.Hmin-V0(2)-Vtend(2))/(p.dt*p.sid);
                QT32 = QV32*T0(3);
                QS32 = QV32*S0(3);
            end
    
        % case where 3 is too thin (and possibly 2, but not 1)
        elseif any(inds==3)
            % box 2 to 3
            QV32 = QV32-(p.W*p.L*p.Hmin-V0(3)-Vtend(3))/(p.dt*p.sid);
            QT32 = QV32*T0(2);
            QS32 = QV32*S0(2);        
             % recalculate tendency for box 2
            Vtend(2) = p.dt*p.sid*(QT0(2)+QVb0(2)+QV32);   
            % if tendency would make box 2 too thin, take from 1
            if V0(2)+Vtend(2)<p.W*p.L*p.Hmin
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