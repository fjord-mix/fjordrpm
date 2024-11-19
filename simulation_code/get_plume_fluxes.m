function [QVp0, QTp0, QSp0, QEp0, QMp0, knb0] = get_plume_fluxes(i, p, s)

% GET_PLUME_FLUXES Compute plume fluxes.
%   [QVp0, QTp0, QSp0, QMp0] = GET_PLUME_FLUXES(i, p, f, s) computes the
%   plume fluxes for the given parameters p, boundary conditions f and
%   solution s at timestep i.

% Get tracer variables at timestep i
H0 = s.H; T0 = s.T(:,i); S0 = s.S(:,i);

% Initialise outputs
[QVp0, QTp0, QSp0, QEp0, QMp0] = deal(zeros(length(p.wp),p.N));
knb0 = zeros(length(p.wp),1);

% Loop over number of plumes
for j = 1:length(p.wp)

    % Get subglacial discharge at timestep i
    Qsg0 = s.Qsg(j,i);
    kgl = s.kgl(j);
    
    if Qsg0~=0 % if there is a plume
        
        % Plume dynamics
        if ~mod(i-1,p.run_plume_every) % if a plume update timestep

            [Qent, Qmelt, knb] = run_plume(j, p, kgl, H0, S0, T0, Qsg0);

        else % otherwise use dynamics from previous time step

            Qent = s.QEp(j,:,i-1);
            Qmelt = s.QMp(j,:,i-1);
            knb = s.knb(j,i-1);

        end
    
        % Compute fluxes in layers from grounding line to neutral buoyancy
        QVp0(j,knb+1:kgl) = -Qent(knb+1:kgl);
        QTp0(j,knb+1:kgl) = QVp0(j,knb+1:kgl).*T0(knb+1:kgl)';
        QSp0(j,knb+1:kgl) = QVp0(j,knb+1:kgl).*S0(knb+1:kgl)';
    
        % Compute fluxes into the neutral buoyancy layer
        Tsg0 = p.l2+p.l3*p.Hgl(j);
        Teff = -p.l/p.cw;
        QVp0(j,knb) = Qsg0 + sum(Qmelt(knb+1:kgl)) - sum(QVp0(j,knb+1:kgl));
        QTp0(j,knb) = Qsg0*Tsg0 + sum(Qmelt(knb+1:kgl))*Teff - sum(QTp0(j,knb+1:kgl));
        QSp0(j,knb) = -sum(QSp0(j,knb+1:kgl));
    
        % Store entrainment, submarine melt flux and neutral buoyancy
        QEp0(j,:) = Qent;
        QMp0(j,:) = Qmelt;
        knb0(j) = knb;
    
    end

end

end