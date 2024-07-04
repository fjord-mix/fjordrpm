function [QpV0, QpT0, QpS0] = get_plume_fluxes(i, p, f, s)

% GET_PLUME_FLUXES Compute plume fluxes.
%   [QpV0, QpT0, QpS0] = GET_PLUME_FLUXES(i, p, f, s) computes the plume
%   fluxes for the given parameters p, boundary conditions f and solution s
%   at timestep i.

% Get tracer variables at timestep i
H0 = s.H; T0 = s.T(:,i); S0 = s.S(:,i);
% Get boundary conditions at timestep i
Qsg0 = f.Qsg(i);

if Qsg0==0 || p.P0==0 
    % If there is no plume, the fluxes are zero by default
    [QpV0, QpT0, QpS0] = deal(0*H0);
else
    % Initialise variables
    [QpV0, QpT0, QpS0] = deal(zeros(p.N, 1));
    
    % Compute the plume properties (flux, salinity, temp) in each zmodel
    % layer and the location of the grounding line and neutral buoyancy
    % box
    [knb, Qp, Sp, Tp] = get_plume_properties(p, s, H0, S0, T0, Qsg0);

    % The flux in boxes below grounding line and above neutral buoyancy are
    % zero, so don't need to be calculated

    % Compute the flux in boxes from the grounding line to below neutral
    % buoyancy
    QpV0(knb+1:s.kgl) = Qp(knb+1:s.kgl)-Qp(knb:s.kgl-1);
    QpT0(knb+1:s.kgl) = QpV0(knb+1:s.kgl).*T0(knb+1:s.kgl);
    QpS0(knb+1:s.kgl) = QpV0(knb+1:s.kgl).*S0(knb+1:s.kgl);

    % Compute the flux into the neutral buoyancy box
    QpV0(knb) = Qp(knb);
    QpT0(knb) = Qp(knb)*Tp(knb);
    QpS0(knb) = Qp(knb)*Sp(knb);
end

end
