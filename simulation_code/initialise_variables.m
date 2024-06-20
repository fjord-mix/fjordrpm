function s = initialise_zmodel_variables(p, a, t)

% INITIALISE_ZMODEL_VARIABLES Initialise zmodel variables with initial
% conditions.
%   S = INITIALISE_ZMODEL_VARIABLES(P, A, T) initialises the zmodel
%   variables using parameters structure P, time T and initial conditions
%   structure A, and returns solution structure S.

% Set the timestep from the inputted time field.
s.dt = t(2:end)-t(1:end-1);

% Set and initialise solution structure fields.
[s.H, s.V, s.T, s.S, s.VT, s.VS,  s.I, ...
    s.QVg, s.QTg, s.QSg, ...
    s.QVs, s.QTs, s.QSs, s.Te, s.Se, s.phi, ...
    s.QVk, s.QTk, s.QSk, ...
    s.QIi, s.QTi, s.QSi, s.QVmi, s.QTmi, s.QSmi, s.M, ...
    s.QVv, s.QTv, s.QSv] = deal(zeros(p.N, length(t)));

% Initialise the layer depths with the given initial conditions.
if p.sill==1
    % If there is a sill, redistribute the layers so that they are roughly
    % the same thickness above and below sill but a box boundary coincides
    % with the sill depth.
    Nabove = round((abs(p.silldepth)/p.H)*p.N);
    Nbelow = p.N-Nabove;
    s.H(:,1) = [(abs(p.silldepth)/Nabove)*ones(1,Nabove),...
        ((p.H-abs(p.silldepth))/Nbelow)*ones(1,Nbelow)];
    % Store the location of the box boundary coinciding with the sill.
    s.ksill = Nabove;
else
    s.H(:,1) = a.H0;
    s.ksill = p.N;
end

% Redistribute the other given initial conditions according to the new
% layer boundaries and then initialise the variables.
ints_old = [0; cumsum(a.H0)];
centres_old = 0.5*(ints_old(1:end-1)+ints_old(2:end));
ints_new = [0;cumsum(s.H(:,1))];
centres_new = 0.5*(ints_new(1:end-1)+ints_new(2:end));

s.T(:,1) = interp1(centres_old,a.T0,centres_new,'linear','extrap');
s.S(:,1) = interp1(centres_old,a.S0,centres_new,'linear','extrap');
s.I(:,1) = interp1(centres_old,a.I0,centres_new,'linear','extrap');

s.V(:,1) = s.H(:,1)'*p.W*p.L; % volume of layers
s.VT(:,1) = s.V(:,1).*s.T(:,1); % heat content
s.VS(:,1) = s.V(:,1).*s.S(:,1); % salt content

end