function [Th,Sh,Vh,VTh,VSh,Hh] = homogenise_layers(V0,T0,S0,VT0, VS0, inds,L,W)

Th = T0;
Sh = S0;
Vh = V0;
VTh=VT0;
VSh=VS0;
% Th(i1) = (T0(i1).*V0(i1) + T0(i2).*V0(i2))./(V0(i1)+V0(i2));
% Sh(i1) = (S0(i1).*V0(i1) + S0(i2).*V0(i2))./(V0(i1)+V0(i2));

Th(inds) = sum(T0(inds).*V0(inds))./sum(V0(inds));
Sh(inds) = sum(S0(inds).*V0(inds))./sum(V0(inds));
Vh(inds) = sum(V0(inds))./length(inds);
VTh(inds) = sum(VT0(inds))./length(inds);
VSh(inds) = sum(VS0(inds))./length(inds);
Hh = Vh./(L*W);

end