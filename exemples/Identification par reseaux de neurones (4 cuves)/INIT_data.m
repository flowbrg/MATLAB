clc; clear; close all;
load("data.mat");
import casadi.*

%% Initialisation des données
% On cherche à ientifier f tq x_k=1 = f(x_k, u_k)

y            = states(:, 2:end)./xmax;                  % x_k+1
x            = [states(:,1:end-1)./xmax; inputs./umax]; % [x_k; u_k]

% Shuffle data
permutation  = randperm(size(x,2));
xperm        = x(:, permutation); yperm= y(:, permutation);

% Separate in data 70/30 proportions
i_test       = floor(0.7*size(x,2));
x_train      = xperm(:, 1:i_test);
y_train      = yperm(:, 1:i_test);
x_valid      = xperm(:,i_test+1:end);
y_valid      = yperm(:,i_test+1:end);