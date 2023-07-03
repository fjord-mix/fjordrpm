%% function to calculate vertical mixing fluxes
function [QVk0,QTk0,QSk0] = mixing_fluxes(H0,T0,S0,QVg0,QVs0,p)

    if p.K0==0

        QVk0 = 0*H0;
        QTk0 = 0*H0;
        QSk0 = 0*H0;

    else
       
        % vector for doing numerical integral
        x = linspace(0,p.L,10);
    
        % loop over interfaces
        for k=1:p.N+p.sill-1
            % velocity in lower box
            ul = (QVg0(k+1)+(QVs0(k+1)-QVg0(k+1))*x/p.L)/(p.W*H0(k+1));
            % velocity in higher box
            uu = (QVg0(k)+(QVs0(k)-QVg0(k))*x/p.L)/(p.W*H0(k));
            % buoyancy jump between boxes
            B = p.g*(p.betaS*(S0(k+1)-S0(k))-p.betaT*(T0(k+1)-T0(k)));

            % an unstable stratification would yield complex velocities if
            % based on the Richardson number
            if B > 0 
                % richardson number
                R = 0.5*(H0(k+1)+H0(k))*B./(ul-uu).^2;
                % vertical entrainment velocity
                % capped at 4e-5
                w = min(p.K0*abs(ul-uu).*R.^(-0.75),p.wmax);
            else
                w = 1.*ones(size(x)).*p.wmax; % so if unstable, we use the min mixing velocity
            end
            % fluxes
            Q(k) = sign(mean(abs(uu)-abs(ul)))*p.W*trapz(x,w);
            QT(k) = (Q(k)>0)*Q(k)*T0(k+1) + (Q(k)<0)*Q(k)*T0(k);
            QS(k) = (Q(k)>0)*Q(k)*S0(k+1) + (Q(k)<0)*Q(k)*S0(k);
        end
    
        % final fluxes
        QVk0 = [Q,0]'-[0,Q]';
        QTk0 = [QT,0]'-[0,QT]';
        QSk0 = [QS,0]'-[0,QS]';
    
    end
       
end