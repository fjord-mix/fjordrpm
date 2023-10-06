function [homogenisation_flag, Ht, Tt, St, Vt, VTt, VSt] = homogenise_thin_layers(Vt, Tt, St, VTt, VSt, Ht, dt, sid, QVg, QVs, QVk, QVb, W, L, N)

homogenisation_flag = false;

% Compute the thickness of each layer at the next timestep.
V_tp1  = Vt+dt*sid*(QVg-QVs+QVk+QVb);
H_tp1 = V_tp1/(W*L);

for k=1:N
    % Only check up to the layer above the bottom.
    if H_tp1(k) <= 0 
        if  k< N &&   H_tp1(k+1) > 0 
            % Else, apply layer homogenisation for each layer k that is a
            % problem.
            homogenisation_flag = true;
            % If layer k+1 is not going to collapse, homogenise k and k+1
            [Tt,St,Vt,VTt,VSt,Ht] = homogenise_layers(Vt, Tt, St, VTt, VSt,[k,k+1],L,W);
        elseif  k < N &&  H_tp1(k+1) < 0 
            % If layer k+1 is also going to collapse, error because we have
            % two adjacent collapsing layers
            error("Error: two adjacent collapsing layers")
        elseif  k == N   && H_tp1(k-1) > 0 
            % in the case that the collapsing layer is adjacent to the
            % bottom/sill, homogenise with the layer above.
            % (This layer is not collapsing or it would have errored already)
            homogenisation_flag = true;
            % If layer k-1 is not going to collapse, homogenise k and k-1
            [Tt,St,Vt,VTt,VSt,Ht] = homogenise_layers(Vt, Tt, St, VTt, VSt,[k,k+1],L,W);
        end
    end
end
end