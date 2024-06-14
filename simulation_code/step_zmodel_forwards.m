function s = step_zmodel_forwards(i, p, s, Q)






    QVg = Q.QVg;
QTg = Q.QTg;
QSg = Q.QSg;

QVs = Q.QVs;
QTs = Q.QTs;
QSs = Q.QSs;

QVk = Q.QVk;
QTk = Q.QTk;
QSk = Q.QSk;

QTi = Q.QTi;
QSi = Q.QSi;
QVmi = Q.QVmi;
QTmi = Q.QTmi;
QSmi = Q.QSmi;

QVv = Q.QVv;
QTv = Q.QTv;
QSv = Q.QSv;

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
