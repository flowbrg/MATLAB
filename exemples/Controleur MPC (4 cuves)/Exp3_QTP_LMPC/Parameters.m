% Author: Yuqi Liu & Lucas Gruss
% Date: 2025-03-03

%% Numerical conditions

clear all;close all; warning off

% set randomization seed
rng(2024-01-25);

% automatic installation of CasADI package
if ~isdir('../casadi-3.6.5')
    %urlwrite('https://github.com/casadi/casadi/releases/download/3.6.5/casadi-3.6.5-windows64-matlab2018b.zip', '../casadi.zip')
    urlwrite('https://github.com/casadi/casadi/releases/download/3.6.5/casadi-3.6.5-linux64-matlab2018b.zip', '../casadi.zip')
    
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


%% MPC/LQR setup

% MPC configuration
MPC = struct();
MPC.dt = 50; % MPC time step
MPC.N = 6; % Steps in MPC prediction horizon
MPC.nu = nu;
MPC.nx = nx;
MPC.xmin = QTP.hmin;
MPC.xmax = QTP.hmax;
MPC.umin = QTP.qmin;
MPC.umax = QTP.qmax;

% Cost function: J = sum(x'Qx + du'Rdu)
MPC.R = eye(2); % Weighting on state vector
MPC.Q = 1e4*eye(4); % Weighting on input vector %
MPC.P = zeros(size(MPC.Q));
MPC.xref = xref; % State reference
MPC.uref = uref; % Input reference

% Case: LMPC (linearization on reference)
% **Define your continous linear model here**
f_lref = casadi.Function('ode',{x, u}, {A*(x-xref)+B*(u-uref)},...
{'state', 'input'}, {'dx'});

MPC.model = integrate_model(f_lref, MPC.dt, 'RK4');

%% Integration
function phi = integrate_model(f, dt, integrator)

    import casadi.*

    names = cellstr(name_in(f));
    x = names{1}; u = names{2};
    state = SX.sym(x, size1_in(f, x), size2_in(f, x));
    controls = SX.sym(u, size1_in(f, u), size2_in(f, u));

    switch integrator
      case 'RK4'
        k1 = f(state, controls);
        k2 = f(state+dt*k1/2.0,controls);
        k3 = f(state+dt*k2/2.0,controls);
        k4 = f(state+dt*k3,controls);
        states_1 = (state+dt*(k1+2*k2+2*k3+k4)/6.0);
        phi = Function('phi', {state,controls},{states_1});

      case 'euler'
        states_1 = (state+dt*f(state, controls));
        phi = Function('phi', {state,controls},{states_1});
    end
end
