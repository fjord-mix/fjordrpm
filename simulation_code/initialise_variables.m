function s = initialise_variables(p, t, f, a)

% INITIALISE_VARIABLES Initialise model variables with initial conditions.
%   s = INITIALISE_VARIABLES(p, t, a) initialises the model variables
%   using parameters structure p, time vector t and initial conditions
%   structure a, and returns solution structure s.

% Set the timestep from the input time field
s.dt = t(2:end)-t(1:end-1);

% Initialise solution structure fields
% Fields with dimensions p.N x 1
[s.H, s.V, s.I] = deal(zeros(p.N, 1));
% Fields with dimensions p.N x length(t)
[s.T, s.S, ...
 s.QVp, s.QTp, s.QSp, s.QMp, ...
 s.QVs, s.QTs, s.QSs, s.Ts, s.Ss, s.phi, ...
 s.QVk, s.QTk, s.QSk, ...
 s.QVi, s.QTi, s.QSi, s.QMi, ...
 s.QVv, s.QTv, s.QSv] = deal(zeros(p.N, length(t)));
% Fields with dimensions 1 x length(t)
[s.Qsg] = deal(zeros(1, length(t)));

% Initialise the layer depths with the given initial conditions
if p.sill==1
    % If there is a sill, redistribute the layers so that they are roughly
    % the same thickness above and below sill but a box boundary coincides
    % with the sill depth
    Nabove = round((p.Hsill/p.H)*p.N);
    Nbelow = p.N-Nabove;
    s.H = [(p.Hsill/Nabove)*ones(Nabove,1);...
           ((p.H-p.Hsill)/Nbelow)*ones(Nbelow,1)];
    % Store the location of the box boundary coinciding with the sill
    s.ksill = Nabove;
else
    s.H = a.H0;
    s.ksill = p.N;
end

% Initialise layer volumes
s.V = s.H*p.W*p.L;

% Get forcings on model layers and at model time steps
[s.Ts, s.Ss, s.Qsg] = bin_forcings(f, s.H, t);

% Find layer with grounding line and store index
ints = cumsum(s.H);
s.kgl = find(ints >= p.Hgl-1e-6, 1);

% Redistribute the other given initial conditions according to the new
% layer boundaries and then initialise the variables
ints_old = [0; cumsum(a.H0)];
centres_old = 0.5*(ints_old(1:end-1)+ints_old(2:end));
ints_new = [0;cumsum(s.H)];
centres_new = 0.5*(ints_new(1:end-1)+ints_new(2:end));

s.T(:,1) = interp1(centres_old,a.T0,centres_new,'linear','extrap');
s.S(:,1) = interp1(centres_old,a.S0,centres_new,'linear','extrap');
s.I = interp1(centres_old,a.I0,centres_new,'linear','extrap');

end