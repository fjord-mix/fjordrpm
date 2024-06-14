function s = step_zmodel_forwards(i, p, s, Q)

    QVg = Q.Qg.V;
QTg = Q.Qg.T;
QSg = Q.Qg.S;

QVs = Q.Qs.V;
QTs = Q.Qs.T;
QSs = Q.Qs.S;

QVk = Q.Qk.V;
QTk = Q.Qk.T;
QSk = Q.Qk.S;

QTi = Q.Qi.T;
QSi = Q.Qi.S;
QVmi = Q.Qi.Vm;
QTmi = Q.Qi.Tm;
QSmi = Q.Qi.Sm;

QVv = Q.Qv.V;
QTv = Q.Qv.T;
QSv = Q.Qv.S;

% Step the temperature, salt, heat content and salt content of the fjord forwards.
    s.V(:,i+1)  = s.V(:,i)+p.dt*p.sid*(QVg-QVs+QVk+QVmi+QVv);
    s.VT(:,i+1) = s.VT(:,i) +p.dt*p.sid*(QTg-QTs+QTk+QTi+QTmi+QTv);
    s.VS(:,i+1) = s.VS(:,i)+p.dt*p.sid*(QSg-QSs+QSk+QSi+QSmi+QSv);
    % compute tracers
    s.H(:,i+1) = s.V(:,i+1)/(p.W*p.L);
    s.T(:,i+1) = s.VT(:,i+1)./s.V(:,i+1);
    s.S(:,i+1) = s.VS(:,i+1)./s.V(:,i+1);


    % Step icebergs forwards.
%     I = I+ (1-p.icestatic)*p.dt*p.sid*((D/(p.W*p.L))*f.xi-M.*I-p.uIce/p.L*I);

end
