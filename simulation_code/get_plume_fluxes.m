function [QVp0, QTp0, QSp0, QMp0] = get_plume_fluxes(i, p, s)

% GET_PLUME_FLUXES Compute plume fluxes.
%   [QVp0, QTp0, QSp0, QMp0] = GET_PLUME_FLUXES(i, p, f, s) computes the
%   plume fluxes for the given parameters p, boundary conditions f and
%   solution s at timestep i.

% Get tracer variables at timestep i
H0 = s.H; T0 = s.T(:,i); S0 = s.S(:,i);
% Get boundary conditions at timestep i
Qsg0 = s.Qsg(i);

if Qsg0==0 % if there is no plume, the fluxes are zero
    
    [QVp0, QTp0, QSp0, QMp0] = deal(0*H0);

else % if Qsg is non-zero, run plume model

    % Initialise variables
    [QVp0, QTp0, QSp0, QMp0] = deal(zeros(p.N, 1));
    
    % Compute the plume properties
    [Qent, Qmelt, knb, Qnb, Snb, Tnb] = run_discreteplume(p, s.kgl, H0, S0, T0, Qsg0);

    % Compute the flux in boxes from the grounding line to below neutral
    % buoyancy
    QVp0(knb+1:s.kgl) = -Qent(knb+1:s.kgl);
    QTp0(knb+1:s.kgl) = QVp0(knb+1:s.kgl).*T0(knb+1:s.kgl);
    QSp0(knb+1:s.kgl) = QVp0(knb+1:s.kgl).*S0(knb+1:s.kgl);

    % Compute the flux into the neutral buoyancy box
    QVp0(knb) = Qnb;
    QTp0(knb) = Qnb*Tnb;
    QSp0(knb) = Qnb*Snb;

    % Store submarine melt flux in solution
    QMp0 = Qmelt;

end

end
