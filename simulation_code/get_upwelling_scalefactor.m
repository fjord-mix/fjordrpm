function scalefac = get_upwelling_scalefactor(p, H0, T0, S0, I0, meltflux, gmelt)

% GET_UPWELLING_SCALEFACTOR Compute the upwelling scalefactor for icebergs
%   scalefac = GET_UPWELLING_SCALEFACTOR(p, H0, T0, S0, I0, meltflux,
%   gmelt) computes the upwelling scalefactor for the amount of iceberg
%   melt mixing between boxes for input parameters p, tracers H0, T0,
%   S0, I0, and iceberg variables meltflux and gmelt.

gk = max(0,[NaN;p.g*(p.betaS*(S0(2:end)-S0(1:end-1))-p.betaT*(T0(2:end)-T0(1:end-1)))]);
lengthfac = (1/p.alphaI^(2/3))*((meltflux./I0).^2./(gmelt.*H0)).^(1/3).*gmelt./gk;
scalefac = 1-exp(-lengthfac);
% Avoid NaN values
scalefac(gk==0) = 1; % overcoming unstable stratification if homogenise_unstable_layers didn't work
scalefac(I0==0) = 0; % no icebergs
% Set no upwelling to atmosphere
scalefac(1) = 0;

end