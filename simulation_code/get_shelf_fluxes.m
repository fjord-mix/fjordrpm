function [QVs0,QTs0,QSs0,Se0,Te0,phi0] = get_shelf_fluxes(H0,T0,S0,zs,Ts,Ss,Qsg0,p)

% GET_SHELF_FLUXES Compute shelf fluxes.
%   [QVS0,QTS0,QSS0,SE0,TE0,PHI0] = GET_SHELF_FLUXES(H0,T0,S0,ZS,TS,SS,QSG0,P)
%   computes the shelf fluxes for the given parameters. 

if p.C0==0
    QVs0 = 0*H0;
    QTs0 = 0*H0;
    QSs0 = 0*H0;
    Se0 = NaN*H0;
    Te0 = NaN*H0;
    phi0 = NaN*H0;

else
    % calculate mean shelf T/S over box model layers above sill
    [gp, phi0, Q] = deal(zeros(1, p.N+p.sill));
    [Te0, Se0, He0] = bin_ocean_profiles_withsill(Ts,Ss,zs,H0',p);

    % compute the location of the below sill, partially occluded and above sill layers
    depthcheck = H0;
    depths = [0; cumsum(H0)];
    for k = 1:p.N + p.sill
        if depths(k) < p.silldepth && depths(k+1) <= p.silldepth
            depthcheck(k) = 1; %'abovesill';
        elseif depths(k) < p.silldepth && depths(k+1) > p.silldepth 
            depthcheck(k) = 2; %'partiallyoccluded';
        elseif depths(k) >= p.silldepth
            depthcheck(k) = 3; %'belowsill';
        end
    end

    % get fjord to shelf reduced gravity
    for k=1:p.N + p.sill
     % no need to split into above, partially occluded and below sill because the shelf
     % conditions have already been computed 
            gp(k) = p.g*(p.betaS*(S0(k)-Se0(k))-p.betaT*(T0(k)-Te0(k)));
    end

    % calculate potentials over above-sill layers
    for k=1:p.N + p.sill
        % first loop is for layer below the surface 
        if depthcheck(k)==1 || depthcheck(k) ==2 % above sill or partially occluded 
            if k==1
                phi0(k) = gp(k)*H0(k)/2;
            else
                phi0(k) = phi0(k-1)+gp(k-1)*H0(k-1)/2+gp(k)*H0(k)/2;
            end
        elseif depthcheck(k) == 3 % below sill
            if depthcheck(k-1) == 2 || depthcheck(k-1)==1 % first below sill layer
                phi0(k) = phi0(k-1)+gp(k-1)*H0(k-1)/2+gp(k)*H0(k);
            else 
                phi0(k) = phi0(k-1)+gp(k)*H0(k)/2;
            end
        end
    end

    % fluxes before barotropic compensation
    for k = 1:p.N + p.sill
            Q(k) = p.C0*p.W*He0(k)*phi0(k)/p.L; % He0 contains information about the gap that the flow can go through 
    end

   % Q = p.C0*p.W*H0(1:p.N).*phi0'/p.L;

    % fluxes after ensuring depth mean = QSg0 when plume is turned on,
    % and 0 when plume is turned off
    if p.P0==0
        Qsg0 = 0;
    end

 
    QVs0 = Q' + He0(1:p.N+p.sill)*(Qsg0-sum(Q))/sum(He0(1:p.N+p.sill));
    % 
    % if p.sill==1
    %     QVs0(p.N+p.sill) = 0;
    %     Se0(p.N+p.sill) = 0;
    %     Te0(p.N+p.sill) = 0;
    % end

    % 
% Check the signs of the below sill fluxes 
% If they disagree with the direction of the flux of the deepest below sill
% box then they become fluxes to the box above/below
% for k = 1:p.N +p.sill
%     if depthcheck(k)==3 % below sill
%         if sign(QVs0(k)) ~= sign(QVs0(end)) 
%             if sign(QVs0(k)) > 0
%                 QVs0(k) = 0;
%                 QVs0()

    % resulting heat/salt fluxes
    QTs0 = (QVs0>0).*QVs0.*T0 + (QVs0<0).*QVs0.*Te0;
    QSs0 = (QVs0>0).*QVs0.*S0 + (QVs0<0).*QVs0.*Se0;
end

end