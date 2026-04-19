% generate data
% Réalise la synthese de commandes MPC (NL et linéarisé) partant du modèle de
% Toricelli.
clear all
%% Préliminaires : Import de Casadi et des fonctions utiles.
close all
addpath("\controllers")
casadi_path = "C:\Users\m25zodro\Documents\TD_MCA_MPC2student\casadi-3.7.2";
%casadi_path = "C:\backup\enseignement\TAF ++++\UE E Méthodologie de Commande Avancee\MINIPROJETs 2023\QTP_lq&MPC\MPC_2022to2023_Lucas\casadi-windows-matlabR2016a-v3.5.5"
if exist(casadi_path) ~= 7
    error('Indiquez un chemin valide vers votre installation de Casadi. <a href="https://web.casadi.org/get/">Instructions pour installation</a>.');
end
addpath(casadi_path);
import casadi.*

%% Synthèse de contrôleurs MPC
%%% Modélisation du processus (modèle interne du contrôleur MPC)
% Le contrôleur MPC utilise un modèle du processus afin de formaliser un
% problème d'optimisation. Les fonctions suivantes nous sont utiles :
% * find_equilibrium: permet de trouver l'équilibre (états et commandes) pour
% des références de hauteur 1 et 2 données
% * QTP_toricelli: renvoie une fonction Casadi qui modélise le système 4
% réservoirs à temps continu. En option, on peut préciser les paramètres du
% modèle si on veut choisir des valeurs différentes des valeurs par défaut.
% * discretize_model: discrétise le modèle avec une méthode de Runge-Kutta et
% un pas d'échantillonage donné. Ce pas correspond à la fréquence de
% fonctionnement de la MPC.

dt = 1;    % temps d'échantillonage (de 1 à 5 sec ?)
f = QTP_toricelli();                  % modèle à temps continu
phi = discretize_model(f, 'dt', dt);   % modèle à temps discret pour des commandes de 5s

% contraintes sur l'état et sur les commandes
xmin = [0 0 0 0]';
xmax = [1.36 1.36 1.3 1.3]';
umin = [0 0]';
umax = [3.26 4]'/3600;

%%% Définition des paramètres communs aux différentes implémentations
N = floor(5/dt); % horizon de commande : Nxdt de l'ordre de 50s, soit N=50 pour dt=1 et N=10 pour dT=5
Q = eye(4); % pondération du cout sur l'écart à la référence
R = 10000000*eye(2); % pondération du cout d'effort de commande
                 % => numeriquement, on a A=phi.jacobian(xref); mais ici, le linéarisé est
                 % codé en dur ci-dessous (pour Xref=findequilibrium(|0,6 ; 0,8])
                 % permet ici de calculer le cout terminal de la MPC avec Lyapunov discret
PHI = phi.jacobian();

%% Generate equilibrium points
x12ref = [0.6 0.8;
          0.8 0.6;
          %0.3 0.9;
          %0.9 0.5;
          %1.0 0.4;
          %0.3 1.1;
          %1 1.15;
          0.43 0.7;
          %0.2 0.7;
          0.7 0.7;
          0.4 0.4;
          0.55 0.2;
          1.2 1.2;
          0.8 0.8;
          0.7 0.9;];

x12ref = 0.95*xmax(1:2)'.*rand(100, 2);

% get full state and control inputs in steady-state
Xref = [];
Uref = [];
for k=1:size(x12ref, 1)
    [temp_xref, temp_uref] = find_equilibrium(x12ref(k, :));

    % filter out impossible references
    if any([(temp_uref > umax)' (temp_uref < umin)' (temp_xref > xmax)' (temp_xref<xmin)'])
        continue;
    end

    Xref = [Xref temp_xref];
    Uref = [Uref temp_uref];
end


Tf = 400;
x0 = [0.65 0.65 0.65 0.65]'; % état initial
t = 0:1:Tf;
%X = [x0];
%U = [];
for k=27:size(Xref, 2)
    phi_eq = full(PHI(Xref(:, k), Uref(:, k), []));
    B = phi_eq(1:4, 5:6);
    A = phi_eq(1:4, 1:4);
    P = dlyap(A, Q);     % cout terminal de la MPC avec Lyapunov discret

    try
        nmpc = MPC(phi, N, Q, R,P,xmin,xmax,umin,umax); % nonlinear MPC
        [U_nmpc, X_nmpc, time_nmpc, pwin] = run_QTP_with_controller(@nmpc.solve, X(:, end), Xref(:, k), t, dt, []);
        close(pwin)
        U = [U U_nmpc];
        X = [X X_nmpc];
    catch
        warning('Could not solve MPC')
    end
end
