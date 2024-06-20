function [QpV0, QpT0, QpS0] = get_plume_fluxes(i, p, f, s)

% GET_PLUME_FLUXES Compute plume fluxes for the zmodel.
%   [QPV0, QPT0, QPS0] = GET_PLUME_FLUXES(I, P, F, S) computes the plume
%   fluxes for the given parameters P, boundary conditions F and solution S
%   at time I.

% Get tracer variables at timestep i.
H0 = s.H(:,i); T0 = s.T(:,i); S0 = s.S(:,i);
% Get boundary conditions at timestep i.
Qsg0 = f.Qsg(i);

if Qsg0==0 || p.P0==0 
    % If there is no plume, the fluxes are zero by default.
    [QpV0, QpT0, QpS0] = deal(0*H0);
else
    % Initialise variables.
    [QpV0, QpT0, QpS0] = deal(zeros(p.N, 1));
    
    % Compute the plume properties (flux, salinity, temp) in each zmodel
    % layer and the location of the grounding line and neutral buoyancy
    % box.
    [kgl, knb, Qp, Sp, Tp] = get_plume_properties(p, H0, S0, T0, Qsg0);

    % The flux in boxes below grounding line and above neutral buoyancy are
    % zero, so don't need to be calculated.

    % Compute the flux in boxes from the grounding line to below neutral
    % buoyancy.
    QpV0(knb+1:kgl) = Qp(knb+1:kgl)-Qp(knb:kgl-1);
    QpT0(knb+1:kgl) = QpV0(knb+1:kgl).*T0(knb+1:kgl);
    QpS0(knb+1:kgl) = QpV0(knb+1:kgl).*S0(knb+1:kgl);

    % Compute the flux into the neutral buoyancy box.
    QpV0(knb) = Qp(knb);
    QpT0(knb) = Qp(knb)*Tp(knb);
    QpS0(knb) = Qp(knb)*Sp(knb);
end


end
