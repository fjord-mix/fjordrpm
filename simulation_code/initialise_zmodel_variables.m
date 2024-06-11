function [H, V, T, S, I, VT, VS,...
    QVg, QTg, QSg, QVs, QTs, QSs,...
    QVk, QTk, QSk, QVmi, QTmi, QSmi,QIi,QTi,QSi,...
    QVv,QTv,QSv,Te,Se, phi, M, p] = initialise_zmodel_variables(p, f, a, t)

[H, V, T, S, I, M] = deal(zeros(p.N, length(t)));

[QVg,QTg,QSg,...
    QVs,QTs,QSs,...
    QVk,QTk,QSk,...
    QVmi,QTmi,QSmi,...
    QIi,QTi,QSi,...
    QVv,QTv,QSv,...
    Te,Se, phi] = deal(zeros(p.N,length(t)-1));

% [I, M] = deal(zeros(length(f.zi),length(t)-1));

%% Initialise variables according with the boundary and initial conditions.
H(:,1) = a.H0;
if p.sill==1
    if p.fixedthickness==1
        % If the layers are fixed thickness, distribute layers so that
        % they are roughly the same thickness above and below sill
        % but a box boundary coincides with the sill depth
        Nabove = round((abs(p.silldepth)/p.H)*p.N);
        Nbelow = p.N-Nabove;
        H(:,1) = [(abs(p.silldepth)/Nabove)*ones(1,Nabove),...
                ((p.H-abs(p.silldepth))/Nbelow)*ones(1,Nbelow)];
        p.ksill = Nabove;
    else 
        % If the layers are variable thickness, add an extra layer for the
        % sill.
       H(:,1) = [(abs(p.silldepth)/p.N)*ones(1,p.N),p.H-abs(p.silldepth)];
    end
else
    p.ksill = p.N;
end


V(:,1) = H(:,1)'*p.W*p.L; % volume of layers
T(:,1) = a.T0; % temperature
S(:,1) = a.S0; % salinity
I(:,1) = a.I0; % iceberg concentration
if p.icestatic
    % a.I0 is the surface area of icebergs in a box
    % so either use the idealised expression here or load from file
    I(:,1) = p.A0*p.if(p.nu0, abs(p.zgl), -cumsum(H(:,1))+H(:,1)/2);    
end
VT(:,1) = V(:,1).*T(:,1); % heat content
VS(:,1) = V(:,1).*S(:,1); % salt content

end