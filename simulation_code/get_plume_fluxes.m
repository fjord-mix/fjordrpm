function [QVp0, QTp0, QSp0, QMp0] = get_plume_fluxes(i, p, s)

% GET_PLUME_FLUXES Compute plume fluxes.
%   [QVp0, QTp0, QSp0, QMp0] = GET_PLUME_FLUXES(i, p, f, s) computes the
%   plume fluxes for the given parameters p, boundary conditions f and
%   solution s at timestep i.

% Get tracer variables at timestep i
H0 = s.H; T0 = s.T(:,i); S0 = s.S(:,i);

% Initialise outputs
[QVp0, QTp0, QSp0, QMp0] = deal(zeros(length(p.wp),p.N));

% Loop over number of plumes
for j = 1:length(p.wp)

    % Get subglacial discharge at timestep i
    Qsg0 = s.Qsg(j,i);
    kgl = s.kgl(j);
    
    if Qsg0~=0 % if there is a plume
        
        % Compute the plume properties
        [Qent, Qmelt, knb, Qnb, Snb, Tnb] = run_discreteplume(j, p, kgl, H0, S0, T0, Qsg0);
    
        % Compute the flux in boxes from the grounding line to below neutral
        % buoyancy
        QVp0(j,knb+1:kgl) = -Qent(knb+1:kgl);
        QTp0(j,knb+1:kgl) = QVp0(j,knb+1:kgl).*T0(knb+1:kgl)';
        QSp0(j,knb+1:kgl) = QVp0(j,knb+1:kgl).*S0(knb+1:kgl)';
    
        % Compute the flux into the neutral buoyancy box
        QVp0(j,knb) = Qnb;
        QTp0(j,knb) = Qnb*Tnb;
        QSp0(j,knb) = Qnb*Snb;
    
        % Store submarine melt flux in solution
        QMp0(j,:) = Qmelt;
    
    end

end

end