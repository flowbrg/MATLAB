% Paramètres
M = 1.0; m = 0.1; l = 0.5; g = 9.81;
dt = 0.02; % pas d'échantillonnage [s]

% Matrices du système linéarisé
A = [0, 1, 0, 0;
    0, 0, -m*g/M, 0;
    0, 0, 0, 1;
    0, 0, (M+m)*g/(M*l), 0];
B = [0; 1/M; 0; -1/(M*l)];

% Condition initiale
theta0 = 10;            % [degrés]
theta0 = pi*theta0/180; % [rad]