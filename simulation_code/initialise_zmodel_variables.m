function s = initialise_zmodel_variables(p, f, a, t)

% initialise variable fieldnames
vars = {'H' 'V' 'T', 'S', 'QVg', 'QTg', 'QSg', 'QVk',  'QVk', 'QTk', 'QSk',...
    'QVmi', 'QTmi', 'QSmi',...
    'QIi', 'QTi', 'QSi',...
    'QVv', 'QTv', 'QSv',...
    'Te', 'Se', 'phi', 'I', 'M', 'VT', 'VS'};

for index = 1 : length(vars)
    s.(vars{index}) = index;  
end

[s.H, s.V, s.T, s.S, s.VT, s.VS] = deal(zeros(p.N, length(t)));

[s.QVg, s.QTg, s.QSg,...
    s.QVs, s.QTs, s.QSs,...
    s.QVk, s.QTk, s.QSk,...
    s.QVmi, s.QTmi, s.QSmi,...
    s.QIi, s.QTi, s.QSi,...
    s.QVv, s.QTv, s.QSv,...
    s.Te, s.Se, s.phi, s.I, s.M] = deal(zeros(p.N, length(t)-1));

% [s.I, s.M] = deal(zeros(length(f.zi), length(t)-1));


% [I, M] = deal(zeros(length(f.zi),length(t)-1));

%% Initialise variables according with the boundary and initial conditions.
s.H(:,1) = a.H0;
%s.V(:,1) = a.H0'*p.W*p.L; % volume of layers
s.T(:,1) = a.T0; % temperature
s.S(:,1) = a.S0; % salinity
s.I(:,1) = a.I0; % iceberg concentration

s.ksill = p.N;
if p.sill==1
    if p.fixedthickness==1
        % If the layers are fixed thickness, redistribute layers so that
        % they are roughly the same thickness above and below sill
        % but a box boundary coincides with the sill depth
        Nabove = round((abs(p.silldepth)/p.H)*p.N);
        Nbelow = p.N-Nabove;
        s.H(:,1) = [(abs(p.silldepth)/Nabove)*ones(1,Nabove),...
                ((p.H-abs(p.silldepth))/Nbelow)*ones(1,Nbelow)];
        s.ksill = Nabove;
    else 
        % If the layers are variable thickness, add an extra layer for the
        % sill.
       s.H(:,1) = [(abs(p.silldepth)/p.N)*ones(1,p.N),p.H-abs(p.silldepth)];
    end
else
    s.ksill = p.N;
end


s.V(:,1) = s.H(:,1)'*p.W*p.L; % volume of layers
%T(:,1) = a.T0; % temperature
%S(:,1) = a.S0; % salinity
%I(:,1) = a.I0; % iceberg concentration

% If the layers are fixed thickness, redistribute initial
% conditions to be the same as the box thicknesses
if p.sill ==1 
    if p.fixedthickness == 1
        % redistribute the layers 
        s.T(:,1) = s.T(:,1).*s.H(:,1)./a.H0;
        s.S(:,1) = s.S(:,1).*s.H(:,1)./a.H0;
        s.I(:,1) = s.I(:,1).*s.H(:,1)./a.H0;
    end
end       
s.VT(:,1) = s.V(:,1).*s.T(:,1); % heat content
s.VS(:,1) = s.V(:,1).*s.S(:,1); % salt content

end