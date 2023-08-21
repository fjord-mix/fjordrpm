function [QVb0, QTb0, QSb0] = get_hmin_fluxes(QV0,QVb0,QTb0,QSb0,V0,H0,T0,S0,p)

    % update the total fluxes vector for possible nudging
    QV0 = QV0 + QVb0;

    % routine to avoid going below minimum thickness
    Vtend = p.dt*p.sid*QV0;
    mininds = find(V0+Vtend<p.W*p.L*p.Hmin);    
    QVmin = zeros(p.N+p.sill,1);            

    % for each box about to go below min thickness, find the nearest
    % (and thickest) box that is not about to go below min thickness 
    % and take volume from that    
    iter_count=1; % we add this counter to prevent an infinite loop
    while ~isempty(mininds) && iter_count < 100
        Has = H0(1:p.N); % thicknesses of layers above the sill
        i_from_box = zeros(size(QVmin)); % we will need to keep track of which box we take the water from for nudging
        for i=1:length(mininds)
        % for i=length(mininds):-1:1
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
        iter_count = iter_count+1;
    end
    if iter_count==100
        disp('did not nudge')
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