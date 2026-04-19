clc; close all; clear;

K = 1;
w_n = 1;

% Temps d'échantillonage
Ts = 0.01;
Np = 1e-10;

% Initialisation
ksi_nom = .2;
x_init = [0;0;ksi_nom];

A = [0 1; -w_n^2 -2*ksi_nom*w_n];
B = [0; K*w_n^2];
C = [1 0];

% Pondérations
Tf = .001;
Rf = 1;
Qf = blkdiag(Tf*eye(2),10);