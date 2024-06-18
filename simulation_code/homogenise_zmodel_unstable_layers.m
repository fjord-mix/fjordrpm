function s = homogenise_zmodel_unstable_layers(i, p, s)

% HOMOGENISE_ZMODEL_UNSTABLE_LAYERS Homogenise two unstably stratified
% layers.
%   S = HOMOGENISE_ZMODEL_UNSTABLE_LAYERS(I, P, S) takes input timestep I,
%   parameters structure P, and solution structure S and computes the
%   buoyancy jump between layers in the solution for that timestep,
%   homogenises any layers that are unstably stratified and returns an
%   updated solution structure S.

% Initialise variables that will be used in a loop- faster than looking
% through the structure in a loop each time.
 V = s.V(:,i);
 T = s.T(:,i);
 S = s.S(:,i);
 VT = s.VT(:,i);
 VS = s.VS(:,i);

% Conmpute the buoyancy jump between boxes.
B = p.g*(p.betaS*(S(2:end)-S(1:end-1))-p.betaT*(T(2:end)-T(1:end-1)));

if any(B < 0)
    % Find the indicies of negative buoyancy entries.
    inx = find(B < 0); 
    for k = 1:length(inx)
        inds = [inx(k), inx(k)+1];
        T(inds) = sum(T(inds).*V(inds))./sum(V(inds));
        S(inds) = sum(S(inds).*V(inds))./sum(V(inds));
        % Recompute heat and salt content.
        VT(inds) = V(inds).*T(inds);
        VS(inds) = V(inds).*S(inds);
    end
end

end

