function [QVp0, QTp0, QSp0, QEp0, QMp0, knb0] = get_plume_fluxes(i, p, s)

% GET_PLUME_FLUXES Compute plume fluxes.
%   [QVp0, QTp0, QSp0, QEp0, QMp0, knb0] = GET_PLUME_FLUXES(i, p, s)
%   computes the plume fluxes for the given parameters p and solution s at
%   timestep i.

% get tracer variables at timestep i
H0 = s.H; T0 = s.T(:,i); S0 = s.S(:,i);

% initialise outputs
[QVp0, QTp0, QSp0, QEp0, QMp0] = deal(zeros(length(p.Wp),p.N));
knb0 = zeros(length(p.Wp),1);

% loop over number of plumes
for j = 1:length(p.Wp)

    % get subglacial discharge at timestep i
    Qsg0 = s.Qsg(j,i);
    kgl = s.kgl(j);
    
    if Qsg0~=0 % if there is a plume
        
        % plume dynamics
        if ~mod(i-1,p.run_plume_every) || s.Qsg(j,i-1)==0 % if a plume update timestep

            [Qent, Qmelt, knb] = run_plume(j, p, kgl, H0, S0, T0, Qsg0);

        else % otherwise use dynamics from previous time step

            Qent = s.QEp(j,:,i-1);
            Qmelt = s.QMp(j,:,i-1);
            knb = s.knb(j,i-1);

        end
    
        % compute fluxes in layers from grounding line to neutral buoyancy
        QVp0(j,knb+1:kgl) = -Qent(knb+1:kgl);
        QTp0(j,knb+1:kgl) = QVp0(j,knb+1:kgl).*T0(knb+1:kgl)';
        QSp0(j,knb+1:kgl) = QVp0(j,knb+1:kgl).*S0(knb+1:kgl)';
    
        % compute fluxes into the neutral buoyancy layer
        Tsg0 = p.l2+p.l3*p.Hgl(j);
        Teff = -p.l/p.cw;
        QVp0(j,knb) = Qsg0 + sum(Qmelt(knb+1:kgl)) - sum(QVp0(j,knb+1:kgl));
        QTp0(j,knb) = Qsg0*Tsg0 + sum(Qmelt(knb+1:kgl))*Teff - sum(QTp0(j,knb+1:kgl));
        QSp0(j,knb) = -sum(QSp0(j,knb+1:kgl));
    
        % store entrainment, submarine melt flux and neutral buoyancy
        QEp0(j,:) = Qent;
        QMp0(j,:) = Qmelt;
        knb0(j) = knb;
    
    end

end

end