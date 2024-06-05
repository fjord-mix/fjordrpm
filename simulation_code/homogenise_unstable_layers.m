function [T, S, VT, VS ] = homogenise_unstable_layers(p, V, T, S, VT, VS)

% Conmpute the buoyancy jump between boxes
B = p.g*(p.betaS*(S(2:end)-S(1:end-1)) - p.betaT*(T(2:end)-T(1:end-1)));

if any(B < 0)
    inx = find(B < 0); % indicies of negative buoyancy entries
    for k = 1:length(inx)
        inds = [inx(k), inx(k)+1];
        T(inds) = sum(T(inds).*V(inds))./sum(V(inds));
        S(inds) = sum(S(inds).*V(inds))./sum(V(inds));
        % Recompute heat and salt content
        VT(inds) = V(inds).*T(inds);
        VS(inds) = V(inds).*S(inds);
    end
end

end

