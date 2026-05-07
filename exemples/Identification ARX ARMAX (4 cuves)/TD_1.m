clear; close all; clc;

%% Paramètres
N = 2047;
Ts = 1;

%% Système réel (ordre 2)
A = [1 -1.5 0.7];
B = [0 0.5 0.2];    % Delay 1
G = tf(B, A, Ts);

%% Signal PRBS
u = idinput(N, 'prbs');

%% Bruit sur l'entrée
sigma_u = 0.1;
u_bruite = u +sigma_u*randn(N,1);

%% Simulation sur la sortie
y_sans_bruit = filter(B,A,u_bruite);

%% Bruit sur la sortie
sigma_y = 0.1;
y = y_sans_bruit + sigma_y*randn(N,1);

%% Données d'identification
% IMPORTANT: on identifie avec l'entrée réellement appliquée
data = iddata(y, u_bruite, Ts);

%% Estimation ARX
na = 2; nb = 2; nk = 1;
modele_arx = arx(data, [na nb nk], 'focus', 'simulation'); %'prediction'

%% Comparaison
compare(data, modele_arx);

%% Affichage des paramètres
disp('Paramètres estimés :')
modele_arx

%% prediction 1-step ahead
y_pred = predict(modele_arx,data, 1);
y_hat_pred = y_pred.OutputData;

% Simuation libre
y_sim = sim(modele_arx, u_bruite); % Entrée utilisée

%% Calcul des indicateurs
% e = y - y_hat_pred;
e = y - y_sim;

MSE = mean(e^2);
RMSE = sqrt(MSE);

NRMSE_std = RMSE / std(y);
NRMSE_range = RMSE / (max(y)-min(y));


%% Estimation ARX
na = 2; nb = 2; nk = 1;
modele_armax = armax(data, [na nb nk],)





N=length(e);                    % longueur des résidus
maxlag = 50;                    % Nombre de décalages à afficher

% Calcul de l'autocorrélation

% Tracé
figure;
stem(lags, acf, 'filled;')