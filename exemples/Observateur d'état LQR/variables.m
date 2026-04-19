clc; clear; close all;

%==============================
% Initialisation des paramètres
%==============================
m1 = 2.29;
m2 = 2.044;
d1_nom = 3.12;
d2 = 3.73;
k = 400;
ky = 280;
km = 2.96;
% Echantillonage
Ts = .01;
Np = 1e-10;

%================================
% Linéarisation du modèle Simulink
%=================================
sys = linmod('modele_bo');

A=sys.a
B=sys.b
C=sys.c
D=sys.d

%=================
% Retour d'état LQ
%=================
Tc = .05;
Qc = inv(Tc*gramt(A,B,Tc));
Rc = 1;
K = lqr(A,B,Qc,Rc);
Kr = 1/(C*inv(-A+B*K)*B);
To = Tc/5;
Qo = inv(To*gramt(A',C',To));
Ro = 1;
L = lqr(A',C',Qo,Ro)';
Ki = 0.005;

%============================
% Observateur idéal de Kalman
%============================
Aobs = A-L*C;
Bobs = [B L];
Cobs = eye(4);
Dobs = zeros(4,2);

