function f = QTP_toricelli(options)
%QTP_TORICELLI returns the continuous-time state evolution function F of the
% 4-tanks system, as well as the discretized function PHI, with sampling time DT.
%
% Usage:
%  [F, PHI] = QTP_toricelli(DT)
%  [F, PHI] = QTP_toricelli(DT, NAME, VALUE), where NAME is any of :
%      - 'S' : section of all tanks.
%      - 'g' : gravity (default 9.81)
%      - 'gamma' : valves' position (default [0.3 0.4])
%      - 'a' : section of the discharge of the tanks (default 0.06)

% https://stackoverflow.com/questions/2775263/how-to-deal-with-name-value-pairs-of-function-arguments-in-matlab/60178631#60178631
arguments
    options.S double = 0.06
    options.g double = 9.81
    options.gamma (2, 1) double = [0.3 0.4] % ouverture des vannes
    options.a (4, 1) double = [1.31e-4 1.51e-4 9.27e-5 8.82e-5]  % sections des écoulements
end
    import casadi.*

    % aliases for brevity
    S = options.S;
    g = options.g;
    gamma = options.gamma;
    a = options.a;

    % Model
    state = SX.sym('x', 1, 4); controls = SX.sym('u', 1, 2); rhs = SX.sym('rhs', 1, 4);
    rhs(1) = (-a(1)*sqrt(2*g*state(1))+a(3)*sqrt(2*g*state(3))+gamma(1)*controls(1))/S;
    rhs(2) = (-a(2)*sqrt(2*g*state(2))+a(4)*sqrt(2*g*state(4))+gamma(2)*controls(2))/S;
    rhs(3) = (-a(3)*sqrt(2*g*state(3))+(1-gamma(2))*controls(2))/S;
    rhs(4) = (-a(4)*sqrt(2*g*state(4))+(1-gamma(1))*controls(1))/S;
    f = Function('f', {state, controls}, {rhs}, {'x', 'utest'}, {'xdot'});
end