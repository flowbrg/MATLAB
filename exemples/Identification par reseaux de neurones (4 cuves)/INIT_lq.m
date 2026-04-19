clc; clear; close all;
load("trained_models.mat");
load("data.mat")
import casadi.*

h = 0.01;

function y_hat = forwardStruct(model, x_k, u_k)
    % x_k+1 = x_k + h*g(x_k) + h*B*u_k
    h = 0.01;
    y_hat = x_k + h * forward(model, x_k, u_k); 
end

%% Point de reference
x_ref = [0.6; 0.8; 0.4334; 1.1449];
u_ref = 1e-4*[5.972; 4.505];
input = dlarray([x_ref./xmax; u_ref./umax], 'CB');

%% Calcul des jacobiennes
function [y_hat, jacX, jacU] = jacobian(model, x_k)
    y_hat = forward(model, dlarray(x_k,'CB'));
    jac   = dljacobian(stripdims(y_hat), stripdims(x_k), 1);
    jacX  = jac(:,1:4);
    jacU  = jac(:,5:6);
end

function [y_hat, jacX, jacU] = jacobian_struct(model, x_k, u_k)
    y_hat = forwardStruct(model, x_k, u_k);
    jacX  = dljacobian(stripdims(y_hat), stripdims(x_k), 1);
    jacU  = dljacobian(stripdims(y_hat), stripdims(u_k), 1);
end

% --- Jacobienne ---
[~, dJdx, dJdu] = dlfeval(@jacobian, net_sigmoid_adam, input);
%[~, dJdx, dJdu] = dlfeval(@jacobian_struct, net_struct, input(1:4,:), input(5:6,:));
A = extractdata(dJdx);
B_lin = extractdata(dJdu);
%% Calcul du gain LQ

% Matrices de pondération de Bryson
Q = diag(1 ./ (xmax.^2));  % [4×4]
R = diag(1 ./ (umax.^2));  % [2×2]

% Calcul du gain LQR (système discret)
[K, ~, ~] = dlqr(A, B_lin, Q, R);

abs(eig(A - B_lin * K)) % <1

%% Simulation en boucle sur le modèle réseau de neurones

N      = 3000;                % nombre de pas
x_sim  = zeros(4, N+1);
u_sim  = zeros(2, N);
x_sim(:,1) = 0.65*ones(4,1);             % condition initiale

for k = 1:N
    % Loi de commande LQR
    u_sim(:,k) = (u_ref./umax - K * (x_sim(:,k)./xmax- x_ref./xmax)).*umax;

    % Saturation optionnelle
    u_sim(:,k) = max([0; 0], min(umax, u_sim(:,k)));

    % Fonction casadi
    x_next = phi(x_sim(:,k), u_sim(:,k));
    x_sim(:,k+1) = full(x_next);
end

% Affichage
t = (0:N) * h;
figure
subplot(2,1,1)
plot(t, x_sim')
legend('x1','x2','x3','x4')
xlabel('Temps (s)'); ylabel('États')
title('États en boucle fermée (casadi)')
yline(x_ref', '--', {'x1_{ref}','x2_{ref}','x3_{ref}','x4_{ref}'})

subplot(2,1,2)
plot(t(1:end-1), u_sim')
legend('u1','u2')
xlabel('Temps (s)'); ylabel('Commandes')
title('Commandes LQR')