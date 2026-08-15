%% Visualisation de la zone — Bateau 2D
clear; clc; close all;

% --- Paramètres de la scène ---
R_obs = 0.2;    % rayon des poteaux
ecart = 0.3;    % distance de sécurité aux obstacles

% Départ et cible
start_pos  = [0, 0];
target_pos = [10, 10];

% Centres des obstacles
obs = [4.0, 8.0; 6.0, 10; 5.0, 0.0; 5.0, 0.8; 5.0, 1.6; 5.0, 2.4;
    5.0, 3.2; 5.0, 4.0; 8.0, 8.0; 6.0, 5.0; 7.0, 9.0; 6.0, 6.0];
n_obs = size(obs, 1);

% Paramètres physiques du bateau - valeurs arbitraires
m = 0.5;     % Masse
f = 0.2;    % Coefficient de frottement
I = 1;      % Inertie de lacet
G = 0.2;    % Constante de couple de la gouverne, fonction de la surface, longueur, etc