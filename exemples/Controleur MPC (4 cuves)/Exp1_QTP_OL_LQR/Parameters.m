% Author: Yuqi Liu & Lucas Gruss
% Date: 2025-03-03

%% Numerical conditions

clear all;close all; warning off

% set randomization seed
rng(2024-01-25);

% automatic installation of CasADI package
if ~isdir('../casadi-3.6.5')
    urlwrite('https://github.com/casadi/casadi/releases/download/3.6.5/casadi-3.6.5-windows64-matlab2018b.zip', '../casadi.zip')
    %urlwrite('https://github.com/casadi/casadi/releases/download/3.6.5/casadi-3.6.5-linux64-matlab2018b.zip', '../casadi.zip')
    
    unzip('../casadi.zip', '../casadi-3.6.5')
end
addpath('../casadi-3.6.5')

%% QTP system setup

% System parameters
QTP = {};
QTP.S = 0.06;        % tanks section
QTP.g = 9.81;        % gravity
QTP.gamma = [0.3, 0.4]';  % position of the valves
QTP.a = [1.31e-4, 1.51e-4, 9.27e-5, 8.82e-5]';  % sections of escapments
QTP.ode = @(x, u) (-QTP.a.*sqrt(2*QTP.g*x(1:4)) ...
                           + [QTP.a(3:4).*sqrt(2*QTP.g*x(3:4)); 0; 0] ...
                           + [QTP.gamma(1);QTP.gamma(2); 1-QTP.gamma(2);1-QTP.gamma(1)]/3600.* ...
                           [u(1); u(2); u(2); u(1)])/QTP.S; % ODE of system dynamic

% Simulation setup
QTP.h0 = [0.001, 1.2, 0.65, 0.65]'; % Initial condition
QTP.hmax = [1.36, 1.36, 1.36, 1.36]'; % State/output upper limit
QTP.hmin = [0, 0, 0, 0]'; % State/output lower limit
QTP.qmax = [3.26, 4]'; % Input upper limit
QTP.qmin = [0, 0]'; % Input lower limit

%% Find root for reference signal

nx = 4; % Number of state
nu = 2; % Number of input
href = [0.8; 0.6]; % Choice of reference for heights 1 & 2.

% Define system symbolically for rootfinder
x = casadi.MX.sym('state', nx, 1); % Define symbol for state
u = casadi.MX.sym('input', nu, 1); % Define symbol for input
f = casadi.Function('ode', {x, u}, {QTP.ode(x, u)}, {'state', 'input'}, {'dx'}); % Define system dynamic

% Find steady state : (xref, uref) such that dx = 0, i.e f(xref, uref) = 0.
xueq = casadi.MX.sym('xueq', 4, 1); % Define decision variable of optimization
rfp = struct('x', xueq, 'g', f([href;xueq(1:2)], xueq(3:4))); % Define optimization problem
s = casadi.rootfinder('solver', 'nlpsol', rfp, struct('nlpsol', 'ipopt')); % Define nlp solver
sol = full(s([href; zeros(2, 1)], [])); % Solve nlp with arg1 - initial guess on xueq:[0.8, 0.6, 0, 0]; arg2 - parameter vector

% Root result
xref = [href; sol(1:2)]; % Stable state reference
uref = sol(3:4); % Stable input reference

% Linearization on equilibrium point
PHI = f.jacobian(); % Define function return gradient of system
[A, B] = PHI(xref, uref, []); % Define problem: get Jacobien on steady state
A = full(A);
B = full(B);

% Alternative rootfinder
% M = [QTP.a(3)*sqrt(2*QTP.g)/QTP.S,0,QTP.gamma(1)/QTP.S,0;
%      0,QTP.a(4)*sqrt(2*QTP.g)/QTP.S,0,QTP.gamma(2)/QTP.S;
%      -QTP.a(3)*sqrt(2*QTP.g)/QTP.S,0,0,(1-QTP.gamma(2))/QTP.S;
%      0,-QTP.a(4)*sqrt(2*QTP.g)/QTP.S,(1-QTP.gamma(1))/QTP.S,0];
% N = [QTP.a(1)*sqrt(2*QTP.g*href(1))/QTP.S;
%      QTP.a(2)*sqrt(2*QTP.g*href(2))/QTP.S;
%      0;
%      0;];
% eq = inv(M)*N;
% xref = [href;eq(1:2).^2];
% uref = eq(3:4);


%% LQR setup
% ================================================
% Case LQR
% ================================================

% Initialisation des paramètres
a1 = 1.31e-4;
a2 = 1.51e-4;
a3 = 9.27e-5;
a4 = 8.82e-5;
ga = 0.3;
gb = 0.4;
Sc = 0.06;
g = 9.81;

% Points de reference
h1ref = 0.8;
h2ref = 0.6;

qbref = ((1-ga)*a1*sqrt(2*g*h1ref)-ga*a2*sqrt(2*g*h2ref))/(1-ga-gb);
qaref = (-gb*a1*sqrt(2*g*h1ref)+(1-gb)*a2*sqrt(2*g*h2ref))/(1-ga-gb);
h3ref = ((1-gb)*qbref)^2/(2*g*a3^2);
h4ref = ((1-ga)*qaref)^2/(2*g*a4^2);

% Jacobiennes analytiques
Alin = (1/Sc)*sqrt(g/2).*[-a1/sqrt(h1ref) 0 a3/sqrt(h3ref) 0;
    0 -a2/sqrt(h2ref) 0 a4/sqrt(h4ref);
    0 0 -a3/sqrt(h3ref) 0;
    0 0 0 -a4/sqrt(h4ref)];

Blin = (1/(Sc*3600))*[ga 0;
    0 gb;
    0 1-gb;
    1-ga 0];

xmax = 1.36*ones(4,1);
umax = [3.26; 4];

Q = diag(1 ./ (xmax.^2)); 
R = diag(1 ./ (umax.^2));

[K,~,~] = lqr(A, B, Q, R);
