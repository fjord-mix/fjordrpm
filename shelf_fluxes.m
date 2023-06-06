%% function to calculate shelf fluxes
function [QVs0,QTs0,QSs0,Se0,Te0,phi0] = shelf_fluxes(H0,T0,S0,zs,Ts,Ss,Qsg0,p)

    if p.C0==0
        
        QVs0 = 0*H0;
        QTs0 = 0*H0;
        QSs0 = 0*H0;
        Se0 = NaN*H0;
        Te0 = NaN*H0;
        phi0 = NaN*H0;
    
    else
    
        % calculate mean shelf T/S over box model layers
        ints = [0;-cumsum(H0)];
        zs0 = unique(sort([zs;-cumsum(H0)]));        
        Ss0 = interp1(zs,Ss,zs0,'pchip','extrap');
        Ts0 = interp1(zs,Ts,zs0,'pchip','extrap');

        % figure; hold on;
        for k=1:length(ints)-1
            inds = find(zs0<=ints(k) & zs0>=ints(k+1));
            if length(inds) == 1 % if there is only one data point, no need to average it
                Se0(k) = Ss0(inds);
                Te0(k) = Ts0(inds);
            else
                Se0(k) = trapz(zs0(inds),Ss0(inds))/H0(k);
                Te0(k) = trapz(zs0(inds),Ts0(inds))/H0(k);                                
                % plot(Ss0(inds),zs0(inds));
            end
        end
        % box 1 does not quite work with the integral due to really sharp
        % gradients
        Se0(1) = mean(Ss0(zs0<=ints(1) & zs0>=ints(2)));
        Te0(1) = mean(Ts0(zs0<=ints(1) & zs0>=ints(2)));
        % scatter(Se0,-cumsum(H0)); yline(ints,':k','linewidth',0.5);
        % figure; plot(Ts0,zs0); hold on; scatter(Te0,-cumsum(H0)); yline(ints,':k','linewidth',0.5);
        Se0 = Se0';
        Te0 = Te0';        
    
        % get fjord to shelf reduced gravity
        for k=1:length(ints)-1
            gp(k) = p.g*(p.betaS*(S0(k)-Se0(k))-p.betaT*(T0(k)-Te0(k)));
        end
    
        % calculate potentials
        for k=1:3
            if k==1
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