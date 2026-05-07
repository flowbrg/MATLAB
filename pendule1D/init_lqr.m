% Controle 1D d'un pendule en position d'équilibre instable

% Régulateur Gain LQR statique

init_params

% Matrices du système linéarisé
C = [1 0 0 0;   % on mesure x
     0 0 1 0];  % on mesure theta

D = zeros(2,1);

% Critère de Bryson
Q = diag([4, 0.25, 25, 0.25]);
R = 0.0025;

% Gain LQR
K = lqr(A, B, Q, R);

poles_lqr = eig(A - B*K)

% Gain observateur
rank(obsv(A, C))   % =4 donc (A,C) observable

poles_obs = 3 * poles_lqr;

L = place(A', C', poles_obs)';
% L = (lqr(A', C', Qobs, Robs))';  % si on veut une approche LQR double
% eig(A - L*C)

A_obs = A - L*C;
B_obs = [B, L];     % entrée [u; y]
C_obs = eye(4);     % on sort tout l'état estimé
D_obs = zeros(4,3);