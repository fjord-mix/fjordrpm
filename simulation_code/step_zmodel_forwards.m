function Z_ip1 = step_zmodel_forwards(i, p, s, Q_i)

    % Compute the zmodel variables at timestep i+1.
    Z_ip1.V = s.V(:,i) +p.dt*p.sid*(Q_i.QVg-Q_i.QVs+Q_i.QVk+Q_i.QVmi+Q_i.QVv);
    Z_ip1.VT = s.VT(:,i)+p.dt*p.sid*(Q_i.QTg-Q_i.QTs+Q_i.QTk+Q_i.QTi+Q_i.QTmi+Q_i.QTv);
    Z_ip1.VS = s.VS(:,i)+p.dt*p.sid*(Q_i.QSg-Q_i.QSs+Q_i.QSk+Q_i.QSi+Q_i.QSmi+Q_i.QSv);
    Z_ip1.H = Z_ip1.V/(p.W*p.L);
    Z_ip1.T = Z_ip1.VT./Z_ip1.V;
    Z_ip1.S = Z_ip1.VS./Z_ip1.V;

    % Step icebergs forwards.
%     I = I+ (1-p.icestatic)*p.dt*p.sid*((D/(p.W*p.L))*f.xi-M.*I-p.uIce/p.L*I);

end
