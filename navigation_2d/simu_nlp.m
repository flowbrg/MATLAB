%% Optimisation de trajectoire — Bateau 2D avec dérapage
%  CasADi + IPOPT
clear; clc; close all;
import casadi.*

init

%% --- Paramètres du NLP ---
N        = 40;       % nombre d'intervalles
w_ctrl   = 0.01;      % poids régularisation commandes
Tf_init  = 30;        % initialisation Tf (s)
Tf_min   = 5;
Tf_max   = 120;

% Bornes commandes
T_max     = 2.0;
theta_max = pi/3;
u_min     = 1e-2;     % évite singularité alpha
u_max     = 5.0;    % vitesse max theorique = T_max/f 

% Cible
xt = target_pos(1); yt = target_pos(2);

%% --- Déclaration des variables symboliques CasADi ---

% Tf libre
Tf = MX.sym('Tf');
h  = Tf / N;

% États : [x, y, phi, alpha, u, r]  taille 6 x (N+1)
X  = MX.sym('X',  6, N+1);

% Commandes : [T, theta]  taille 2 x N
U  = MX.sym('U',  2, N);

%% --- Fonction de dynamique ---
% Entrée : état x (6x1), commande uc (2x1)
% Sortie : xdot (6x1)

x_s   = MX.sym('x_s',  6);
uc_s  = MX.sym('uc_s', 2);

phi_s   = x_s(3);
alpha_s = x_s(4);
us_s    = x_s(5);
r_s     = x_s(6);
T_s     = uc_s(1);
theta_s = uc_s(2);
beta_s  = phi_s + alpha_s;

f_dyn = Function('f_dyn', {x_s, uc_s}, { ...
    vertcat( ...
        us_s * cos(beta_s), ...
        us_s * sin(beta_s), ...
        r_s, ...
        -(T_s * sin(alpha_s)) / (m * us_s) - r_s, ...
        (T_s * cos(alpha_s) - f * us_s^2) / m, ...
        (G * us_s * sin(alpha_s - theta_s)) / I  ...
    )});

%% --- Construction du NLP ---

J    = 0;           % critère
g    = {};          % contraintes d'égalité (collocation)
g_lb = {};          % bornes inférieures
g_ub = {};          % bornes supérieures

for k = 1:N
    xk  = X(:, k);
    xk1 = X(:, k+1);
    uk  = U(:, k);

    % Critère : distance à la cible + régularisation
    J = J + ((xk(1)-xt)^2 + (xk(2)-yt)^2 + ...
             w_ctrl*(uk(1)^2 + uk(2)^2)) * h;

    % Collocation trapézoïdale
    fk  = f_dyn(xk,  uk);
    fk1 = f_dyn(xk1, uk);
    col = xk1 - xk - (h/2) * (fk + fk1);

    g{end+1}    = col;
    g_lb{end+1} = zeros(6,1);
    g_ub{end+1} = zeros(6,1);
end

% Terme terminal sur la cible
J = J + ((X(1,N+1)-xt)^2 + (X(2,N+1)-yt)^2);

% Contraintes d'obstacles (tous les noeuds)
for k = 1:N+1
    xk = X(:,k);
    for i = 1:n_obs
        dist2 = (xk(1)-obs(i,1))^2 + (xk(2)-obs(i,2))^2;
        g{end+1}    = dist2;
        g_lb{end+1} = (R_obs+ecart)^2;
        g_ub{end+1} = inf;
    end
end

%% --- Conditions aux limites (égalités) ---
% État initial
g{end+1}    = X(:,1) - [0; 0; pi/4; 0; 0.01; 0];
g_lb{end+1} = zeros(6,1);
g_ub{end+1} = zeros(6,1);

% État final : position uniquement
g{end+1}    = X(1:2, N+1) - [xt; yt];
g_lb{end+1} = zeros(2,1);
g_ub{end+1} = zeros(2,1);

%% --- Assemblage vecteur de décision ---
%w     = {Tf, reshape(X, [], 1), reshape(U, [], 1)};
w_vec = vertcat(Tf, vec(X), vec(U)); %vertcat(w{:});
g_vec = vertcat(g{:});
g_lb_vec = vertcat(g_lb{:});
g_ub_vec = vertcat(g_ub{:});

%% --- Bornes sur les variables de décision ---

% Tf
w_lb = Tf_min;
w_ub = Tf_max;

% États X : [x, y, phi, alpha, u, r] x (N+1)
x_lb = [-inf; -inf; -inf; -inf; u_min; -inf];
x_ub = [ inf;  inf;  inf;  inf; u_max;  inf];
w_lb = [w_lb; repmat(x_lb, N+1, 1)];
w_ub = [w_ub; repmat(x_ub, N+1, 1)];

% Commandes U : [T, theta] x N
u_lb = [0;        -theta_max];
u_ub = [T_max;     theta_max];
w_lb = [w_lb; repmat(u_lb, N, 1)];
w_ub = [w_ub; repmat(u_ub, N, 1)];

%% --- Point initial ---
% Interpolation linéaire en position, reste nul
x_init_traj = zeros(6, N+1);
for k = 1:N+1
    s = (k-1)/N;
    x_init_traj(1,k) = s * xt;
    x_init_traj(2,k) = s * yt;
    x_init_traj(3,k) = pi/4;
    x_init_traj(5,k) = 0.5;    % vitesse initiale non nulle
end
u_init = zeros(2, N);
u_init(1,:) = 0.5;             % poussée initiale modérée

w0 = [Tf_init; reshape(x_init_traj, [], 1); reshape(u_init, [], 1)];

%% --- Solveur IPOPT ---
nlp  = struct('x', w_vec, 'f', J, 'g', g_vec);
opts = struct();
opts.ipopt.max_iter        = 2000;
opts.ipopt.tol             = 1e-6;
opts.ipopt.print_level     = 5;

solver = nlpsol('solver', 'ipopt', nlp, opts);

sol = solver('x0',  w0, ...
             'lbx', w_lb, 'ubx', w_ub, ...
             'lbg', g_lb_vec, 'ubg', g_ub_vec);

%% --- Extraction de la solution ---
w_sol  = full(sol.x);

Tf_sol = w_sol(1);
idx_X  = 2 : 6*(N+1)+1;
idx_U  = 6*(N+1)+2 : length(w_sol);

X_sol  = reshape(w_sol(idx_X), 6, N+1);
U_sol  = reshape(w_sol(idx_U), 2, N);
t_sol  = linspace(0, Tf_sol, N+1);

fprintf('Tf optimal : %.2f s\n', Tf_sol);

%% --- Figure 1 : Trajectoire ---
theta_c = linspace(0, 2*pi, 100);
figure(1); hold on; axis equal; grid on;
xlim([-1,12]); ylim([-1,12]);
xlabel('x (m)'); ylabel('y (m)');
title(sprintf('Trajectoire optimale  (Tf = %.1f s)', Tf_sol));

for i = 1:n_obs
    fill(obs(i,1)+R_obs*cos(theta_c), obs(i,2)+R_obs*sin(theta_c), ...
         [0.8 0.2 0.2], 'EdgeColor','k');
end

plot(0, 0, 'gs', 'MarkerSize',12, 'MarkerFaceColor','g', 'DisplayName','Départ');
plot(xt, yt, 'p', 'MarkerSize',16, 'MarkerFaceColor','y', ...
     'MarkerEdgeColor','k', 'DisplayName','Cible');
plot(X_sol(1,:), X_sol(2,:), 'b-', 'LineWidth', 2, 'DisplayName','Trajectoire');

n_arr = 15;
idx   = round(linspace(1, N+1, n_arr));
beta_sol = X_sol(3,:) + X_sol(4,:);
quiver(X_sol(1,idx), X_sol(2,idx), ...
       0.4*cos(beta_sol(idx)), 0.4*sin(beta_sol(idx)), ...
       0, 'k', 'DisplayName','Direction vitesse');
%legend('Location','northwest');

%% --- Figure 2 : États ---
figure(2);
subplot(4,1,1); plot(t_sol, X_sol(5,:), 'b', 'LineWidth',1.5);
ylabel('u (m/s)'); grid on; title('Vitesse longitudinale');

subplot(4,1,2); plot(t_sol, rad2deg(X_sol(3,:)), 'r', 'LineWidth',1.5);
ylabel('\phi (°)'); grid on; title('Cap');

subplot(4,1,3); plot(t_sol, rad2deg(X_sol(4,:)), 'g', 'LineWidth',1.5);
ylabel('\alpha (°)'); grid on; title('Angle de dérapage');

subplot(4,1,4); plot(t_sol, X_sol(6,:), 'm', 'LineWidth',1.5);
ylabel('r (rad/s)'); grid on; title('Vitesse de lacet');
xlabel('t (s)');

%% --- Figure 3 : Commandes ---
t_u = t_sol(1:end-1);
figure(3);
subplot(2,1,1); stairs(t_u, U_sol(1,:), 'b', 'LineWidth',1.5);
ylabel('T (N)'); grid on; title('Poussée');

subplot(2,1,2); stairs(t_u, rad2deg(U_sol(2,:)), 'r', 'LineWidth',1.5);
ylabel('\theta (°)'); grid on; title('Angle de gouverne');
xlabel('t (s)');