function [QVb0,QTb0,QSb0]=get_nudging_fluxes(QVb0,QTb0,QSb0,H0,T0,S0,Ss,zs,p)
    % Preallocate variables
    [zn, QVnudge, QTnudge, QSnudge] = deal(zeros(1, p.N-1));
    % get depth of layer boundaries
    Z0 = cumsum(H0);
    % get desired depth based on salinity boundaries on shelf
    for k = 1:p.N-1
        % desired depth
        [~,minind] = min(abs(Ss-p.Snudge(k)));
        zn(k) = abs(zs(minind));
        % nudging flux
        QVnudge(k) = p.W*p.L*(zn(k)-Z0(k))/(p.trelax*p.sid);
        % associated heat and salt fluxes
        QTnudge(k) = (QVnudge(k)>0)*QVnudge(k)*T0(k+1)+(QVnudge(k)<0)*QVnudge(k)*T0(k);
        QSnudge(k) = (QVnudge(k)>0)*QVnudge(k)*S0(k+1)+(QVnudge(k)<0)*QVnudge(k)*S0(k);
    end

    % update box-to-box vectors
    if p.sill == 1
        QVnudge = [QVnudge,0];
        QTnudge = [QTnudge,0];
        QSnudge = [QSnudge,0];
    end
    QVb0 = QVb0+[QVnudge,0]'-[0,QVnudge]';
    QTb0 = QTb0+[QTnudge,0]'-[0,QTnudge]';
    QSb0 = QSb0+[QSnudge,0]'-[0,QSnudge]';
end