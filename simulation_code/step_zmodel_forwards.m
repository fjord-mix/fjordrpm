function [V, VT, VS, H, T, S, I] = step_zmodel_forwards(p, f, V, VT, VS, I, M, D,...
    QVg, QVs, QVk, QVmi, QVv, ...
    QTg, QTs, QTk, QTi, QTmi, QTv, ...
    QSg, QSs, QSk, QSi, QSmi, QSv)

% Step the temperature, salt, heat content and salt content of the fjord forwards.
    V  = V+p.dt*p.sid*(QVg-QVs+QVk+QVmi+QVv);
    VT = VT +p.dt*p.sid*(QTg-QTs+QTk+QTi+QTmi+QTv);
    VS = VS+p.dt*p.sid*(QSg-QSs+QSk+QSi+QSmi+QSv);
    H = V/(p.W*p.L);
    T = VT./V;
    S = VS./V;


    % Step icebergs forwards.
%     I = I+ (1-p.icestatic)*p.dt*p.sid*((D/(p.W*p.L))*f.xi-M.*I-p.uIce/p.L*I);

end
