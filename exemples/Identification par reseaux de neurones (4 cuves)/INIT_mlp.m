clc; clear; close all;
load("data.mat");
import casadi.*

INIT_data

%% Initialisation du modèle
n_in       = size(x,1);
n_hidden   = 12;
n_out      = size(y,1);
net = dlnetwork;
tempNet = [
    featureInputLayer(n_in,"Name","input")
    fullyConnectedLayer(n_hidden,"Name","hidden_layer")
    fullyConnectedLayer(n_out,"Name","output_layer")];
net_tanh_0 = addLayers(net,[tempNet; tanhLayer("Name","tanh")]);
net_relu_0 = addLayers(net,[tempNet; reluLayer("Name","relu")]);
net_sigmoid_0 = addLayers(net,[tempNet; sigmoidLayer("Name","sigmoid")]);

% clean up helper variable
clear tempNet;

net_relu = initialize(net_relu_0);
net_tanh = initialize(net_tanh_0);
net_sigmoid = initialize(net_sigmoid_0);

%% Fonction d'entrainement gradient

% Loss function
function [loss, gradients] = lossFunction(model, x, y)
    % Forward pass
    y_hat = forward(model, dlarray(x, 'CB'));
    
    % Compute MSE loss
    loss = mse(y_hat, dlarray(y, 'CB'));
    
    % Compute gradients
    gradients = dlgradient(loss, model.Learnables);
end

function net = gradientTraining(net, nb_batch, batch_size, nb_epochs, eta, x, y)
for epoch = 1:nb_epochs
    epoch_loss = 0;
    for i = 1:nb_batch
        idx     = (i-1)*batch_size + 1 : i*batch_size;
        x_batch = dlarray(x(:, idx), 'CB');
        y_batch = dlarray(y(:, idx), 'CB');

        [loss, gradients_relu] = dlfeval(@lossFunction, net, x_batch, y_batch);
        net  = dlupdate(@(w,g) w - eta*g, net, gradients_relu);
        epoch_loss = epoch_loss + extractdata(loss);
    end
    fprintf('Epoch %d/%d - Loss: %.6f\n', epoch, nb_epochs, epoch_loss/nb_batch);
end

end

%% Initialisation de l'entrainement

nb_batch   = 113; % x_train size is 8701 = 7x11x113
batch_size = 77;    % Taille d'un batch
eta        = 1e-2;  % Learning rate
nb_epochs  = 50;    % Nombre d'iterations d'entrainement


% Entrainement
fprintf("\n==== Relu ====\n")
net_relu = gradientTraining(net_relu, nb_batch, batch_size, nb_epochs, eta, x_train, y_train);
fprintf("\n==== Tanh ====\n")
net_tanh = gradientTraining(net_tanh, nb_batch, batch_size, nb_epochs, eta, x_train, y_train);
fprintf("\n=== Sigmoid ===\n")
net_sigmoid = gradientTraining(net_sigmoid, nb_batch, batch_size, nb_epochs, eta, x_train, y_train);


%% Fonction d'entrainement adam

% Adam training
function net = adamTraining(net, nb_batch, batch_size, nb_epochs, eta, x, y)
avg_grad    = [];
avg_sq_grad = [];
for epoch = 1:nb_epochs
    epoch_loss = 0;
    for i = 1:nb_batch
        idx     = (i-1)*batch_size + 1 : i*batch_size;
        x_batch = dlarray(x(:, idx), 'CB');
        y_batch = dlarray(y(:, idx), 'CB');

        [loss, gradients] = dlfeval(@lossFunction, net, x_batch, y_batch);
        [net,avg_grad,avg_sq_grad] = adamupdate(net, gradients, avg_grad, avg_sq_grad, epoch+i-1, eta);
        epoch_loss = epoch_loss + extractdata(loss);
    end

    fprintf('Epoch %d/%d - Loss: %.6f\n', epoch, nb_epochs, epoch_loss/nb_batch);
end
end

%% Entrainement

net_relu_adam = initialize(net_relu_0);
net_tanh_adam = initialize(net_tanh_0);
net_sigmoid_adam = initialize(net_sigmoid_0);

fprintf("\n==== Relu - Adam ====\n")
net_relu_adam = adamTraining(net_relu_adam, nb_batch, batch_size, nb_epochs, eta, x_train, y_train);
fprintf("\n==== Tanh - Adam ====\n")
net_tanh_adam = adamTraining(net_tanh_adam, nb_batch, batch_size, nb_epochs, eta, x_train, y_train);
fprintf("\n=== Sigmoid - Adam ===\n")
net_sigmoid_adam = adamTraining(net_sigmoid_adam, nb_batch, batch_size, nb_epochs, eta, x_train, y_train);


%% Validation du meilleur modèle (net_sigmoid_adam) sur données de validation

function validation(net,x_valid,y_valid)
    % Forward pass on validation data
    x_valid_dl = dlarray(x_valid, 'CB');         % Convert to dlarray (Channel x Batch)
    y_pred_dl  = predict(net, x_valid_dl);
    y_pred     = extractdata(y_pred_dl);          % Convert back to numeric array
    
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
end


fprintf("\n==== Relu - Gradient ====")
validation(net_relu, x_valid, y_valid)
fprintf("\n==== Tanh - Gradient ====")
validation(net_tanh, x_valid, y_valid)
fprintf("\n=== Sigmoid - Gradient ===")
validation(net_sigmoid, x_valid, y_valid)
fprintf("\n==== Relu - Adam ====")
validation(net_relu_adam, x_valid, y_valid)
fprintf("\n==== Tanh - Adam ====")
validation(net_tanh_adam, x_valid, y_valid)
fprintf("\n=== Sigmoid - Adam ===")
validation(net_sigmoid_adam, x_valid, y_valid)

%%
save('trained_models.mat', 'net_relu');
save('trained_models.mat', 'net_tanh', '-append');
save('trained_models.mat', 'net_sigmoid', '-append');
save('trained_models.mat', 'net_relu_adam', '-append');
save('trained_models.mat', 'net_tanh_adam', '-append');
save('trained_models.mat', 'net_sigmoid_adam', '-append');