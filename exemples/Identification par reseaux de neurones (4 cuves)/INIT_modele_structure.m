clc; clear; close all;
load("data.mat");
import casadi.*

INIT_data

%% Initialisation du modèle
n_hidden   = 12;
n_out      = size(y,1);
h = 0.01; % temps d'échantillonage
net        = dlnetwork;

tempNet = [
    featureInputLayer(4,"Name","x_k")
    fullyConnectedLayer(n_hidden,"Name","hidden")
    sigmoidLayer("Name","sigmoid")
    fullyConnectedLayer(n_out,"Name","sortie")];
net = addLayers(net,tempNet);

tempNet = [
    featureInputLayer(2,"Name","u_k")
    fullyConnectedLayer(n_out,"Name","B")];
net = addLayers(net,tempNet);

tempNet = additionLayer(2,"Name","addition");
net = addLayers(net,tempNet);

% clean up helper variable
clear tempNet;

net = connectLayers(net,"sortie","addition/in1");
net = connectLayers(net,"B","addition/in2");
net_struct = initialize(net);

% nouvelle fonction forward
function y_hat = forwardStruct(model, x_k, u_k)
    % x_k+1 = x_k + h*g(x_k) + h*B*u_k
    h = 0.01;
    y_hat = x_k + h * forward(model, x_k, u_k); 
end

%% Entrainement
function [loss, gradients, jac] = lossFunctionStruct(model, x_k, u_k, y)
    % Forward pass
    y_hat = forwardStruct(model, x_k, u_k);
    
    % Compute MSE loss
    loss = mse(y_hat, y);
    
    y_hat = stripdims(y_hat);
    input = stripdims([x_k; u_k]);

    % Compute gradients
    gradients = dlgradient(loss, model.Learnables);

    jac = dljacobian(y_hat, input, 1);
end

% Initialisation
nb_batch   = 113; % x_train size is 8701 = 7x11x113
batch_size = 77;    % Taille d'un batch
eta        = 1e-3;  % Learning rate
nb_epochs  = 20;    % Nombre d'iterations d'entrainement

avg_grad      = [];
avg_sq_grad   = [];
for epoch = 1:nb_epochs
    epoch_loss = 0;
    for i = 1:nb_batch
        idx     = (i-1)*batch_size + 1 : i*batch_size;
        x_batch = dlarray(x_train(1:4, idx), 'CB');  % uniquement x_k
        u_batch = dlarray(x_train(5:6, idx), 'CB');  % uniquement u_k
        y_batch = dlarray(y_train(:, idx), 'CB');

        [loss, gradients, ~] = dlfeval(@lossFunctionStruct, net_struct, x_batch, u_batch, y_batch);
        [net_struct, avg_grad, avg_sq_grad] = adamupdate(net_struct, gradients, avg_grad, avg_sq_grad, epoch+i-1, eta);
        epoch_loss = epoch_loss + extractdata(loss);
    end

    fprintf('Epoch %d/%d — Loss: %.6f\n', epoch, nb_epochs, epoch_loss/nb_batch);
end

%% Validation

% Forward pass on validation data
x_valid_dl = dlarray(x_valid(1:4,:), 'CB');         % Convert to dlarray (Channel x Batch)
u_valid_dl = dlarray(x_valid(5:6,:), 'CB');
y_pred_dl  = forwardStruct(net_struct, x_valid_dl, u_valid_dl);
y_pred = extractdata(y_pred_dl);          % Convert back to numeric array

% Compute residuals
residuals  = y_valid - y_pred;

% --- Metrics ---
% MSE per state
MSE_per_state = mean(residuals.^2, 2);

% RMSE per state
RMSE_per_state = sqrt(MSE_per_state);

% R² (coefficient of determination) per state
SS_res = sum(residuals.^2, 2);
SS_tot = sum((y_valid - mean(y_valid, 2)).^2, 2);
R2_per_state = 1 - SS_res ./ SS_tot;

% Display results
fprintf('\n================== Validation Results  ==================\n');
for i = 1:4
    fprintf('State %d | MSE: %.6f | RMSE: %.6f | R²: %.4f\n', ...
        i, MSE_per_state(i), RMSE_per_state(i), R2_per_state(i));
end
fprintf('=======================================================\n');
fprintf('Mean R²  : %.4f\n', mean(R2_per_state));
fprintf('Mean RMSE: %.6f\n', mean(RMSE_per_state));

save('trained_models.mat', 'net_struct', '-append')
%% Synthèse LQR

% Point de reference
x_ref    = [0.6; 0.8; 0.4334; 1.1449]; % x équilibre
x_ref_dl = dlarray(x_ref./xmax, 'CB'); % x équilibre normalisé
u_ref    = 1e-4*[5.972; 4.505];        % u équilibre
u_ref_dl = dlarray(u_ref./umax, 'CB'); % u équilibre normalisé

% Fonctions nécessaires à dlfeval

function [y_hat, jacX, jacU] = jacobian(model, x_k, u_k)
    y_hat = forwardStruct(model, x_k, u_k);
    jacX  = dljacobian(stripdims(y_hat), stripdims(x_k), 1);
    jacU  = dljacobian(stripdims(y_hat), stripdims(u_k), 1);
end

% Jacobienne par rapport à x
[~, dJdx, dJdu] = dlfeval(@jacobian, net_struct, x_ref_dl, u_ref_dl);
A = extractdata(dJdx);  % [4×4]
B_lin = extractdata(dJdu);  % [4x2]

%% Calcul du gain LQ

% Matrices de pondération de Bryson
% J = int(zQz' + uRu')
Q = diag(1 ./ (xmax.^2));  % [4x4]
R = diag(1 ./ (umax.^2));  % [2x2]

% Calcul du gain LQR (système discret)
[K, S, P] = dlqr(A, B_lin, Q, R);

abs(eig(A - B_lin * K)) % <1

%% Simulation en boucle sur le modèle réseau de neurones

N      = 3000;                % nombre de pas
x_sim  = zeros(4, N+1);
u_sim  = zeros(2, N);
% x_sim(:,1) = x0;            % condition initiale

for k = 1:N
    % Loi de commande LQR
    u_sim(:,k) = (u_ref./umax + K * (x_sim(:,k)./xmax - x_ref./xmax)).*umax;

    % Saturation optionnelle
    u_sim(:,k) = max(umin, min(umax, u_sim(:,k)));

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