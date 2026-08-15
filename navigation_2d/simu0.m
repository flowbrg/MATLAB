%% Simulation navigation 2D basique
% Pas de dérapage
% Le bateau ne se redresse pas quand on place la gouverne en position
% "neutre".
% C'est nul, lancer simu_derapage pour le système "augmenté" de l'angle
% d'incidence.

clear; clc; close all;
init

%% Signaux de commande
% Paramètres des commandes
% Poussée, échelon unitaire
T0  = 1.0; 

% Couple, trapèze
% Montée linéaire de 0 à theta_max entre t1 et t2
% Palier à theta_max entre t2 et t3
% Descente linéaire de theta_max à 0 entre t3 et t4
theta_max   = 3/4;
t1 = 2; t2 = 4; t3 = 6; t4 = 8;

% Définition des commandes
T_cmd       = @(t) T0;  % Echelon unitaire
theta_cmd = @(t) ...
    (t >= t1 & t < t2)  .* (theta_max * (t - t1)/(t2 - t1)) + ...
    (t >= t2 & t < t3)  .* theta_max + ...
    (t >= t3 & t <= t4) .* (theta_max * (t4 - t)/(t4 - t3));

%% Conditions initiales
% [x, y, phi, u, r)
x0 = [0; 0; pi/4; 0; 0];    % Cap initial de 45°, vers la cible

% Intégration ode45, explicit runge-kutta
t_span          = [0,20];

ode_fun         = @(t,x) dynamics(t, x, T_cmd, theta_cmd, m, I, f, G);
options         = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
[t_sol, x_sol]  = ode45(ode_fun, t_span, x0, options);

T_sol           = arrayfun(T_cmd, t_sol);
theta_sol       = arrayfun(theta_cmd, t_sol);

%% Figure 1 : Trajectoire
figure(1); hold on; axis equal; grid on;
xlim([-1, 12]); ylim([-1, 12]);
xlabel('x (m)'); ylabel('y (m)');
plot(x_sol(:,1), x_sol(:,2), 'b-', 'LineWidth', 2, 'DisplayName', 'Trajectoire');

% Cap du bateau
n_arrows = 15;  % Nombre de flèches
idx = round(linspace(1, length(t_sol), n_arrows));
quiver(x_sol(idx,1), x_sol(idx, 2), 0.4*cos(x_sol(idx,3)),  0.4*sin(x_sol(idx,3)), 0, 'k', 'LineWidth', 1.2, 'DisplayName', 'Cap \phi');

legend('Location', 'northwest');

%% Figure 2 : Etats
figure('Name', 'États');

% vitesse longi
subplot(3,1,1);
plot(t_sol, x_sol(:,4), 'b', 'LineWidth', 1.5);
ylabel('v (m/s)'); grid on; title('Vitesse longitudinale');

% cap
subplot(3,1,2);
plot(t_sol, rad2deg(x_sol(:,3)), 'r', 'LineWidth', 1.5);
ylabel('\psi (°)'); grid on; title('Cap');

% vitesse de lacet
subplot(3,1,3);
plot(t_sol, x_sol(:,5), 'm', 'LineWidth', 1.5);
ylabel('r (rad/s)'); grid on; title('Vitesse de lacet');
xlabel('t (s)');

%% Figure 3 : Commandes
figure('Name', 'Commandes');

% Poussée
subplot(2,1,1);
plot(t_sol, T_sol, 'b', 'LineWidth', 1.5);
ylabel('T (N)'); grid on; title('Poussée — échelon');

% Gouverne
subplot(2,1,2);
plot(t_sol, theta_sol, 'r', 'LineWidth', 1.5);
ylabel('u (N.m)'); grid on; title('Couple de gouverne — trapèze');
xlabel('t (s)');

%% Fonction dynamique
function dxdt = dynamics(t, x, T_cmd, theta_cmd, m, I, f, G)
    phi = x(3);
    u   = x(4);
    r   = x(5);

    T       = T_cmd(t);
    theta   = theta_cmd(t);

    dxdt    = zeros(5,1);
    dxdt(1) = u*cos(phi);
    dxdt(2) = u*sin(phi);
    dxdt(3) = r;
    dxdt(4) = (T-f*u)/m;
    dxdt(5) = -(G*u*theta)/I;   % Hyp: theta << 1, sin(theta)=theta
end