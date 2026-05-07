% Controle 1D d'un pendule en position d'équilibre instable

% Régulateur MPC

init_params % paramètres du modèle

init_lqr % initialisation de l'observateur

% Discrétisation exacte (Zero-Order Hold)
sys_c = ss(A, B, eye(4), zeros(4,1));
sys_d = c2d(sys_c, dt, 'zoh');
Ad = sys_d.A;
Bd = sys_d.B;

% Paramètres MPC
N  = 20;           % horizon de prédiction
Q  = diag([4, 0.25, 25, 0.25]);   % on prend meme que LQR
R  = 0.0025;
P  = idare(Ad, Bd, Q, R);         % cout

u_max =  20;   % [N] saturation actionneur
u_min = -20;

% Construction de F et G
nx = 4; nu = 1;
F = zeros(nx*N, nx);
G = zeros(nx*N, nu*N);

Apow = eye(nx);
for i = 1:N
    Apow = Ad * Apow;
    F((i-1)*nx+1 : i*nx, :) = Apow;
    for j = 1:i
        G((i-1)*nx+1 : i*nx, (j-1)*nu+1 : j*nu) = Ad^(i-j) * Bd;
    end
end

% Matrices cout bloc-diagonales
Q_bar = blkdiag(kron(eye(N-1), Q), P);  % P en terminal
R_bar = kron(eye(N), R);

% Pré-calculer le produit constant G'*Q_bar*F
GtQbarF = G' * Q_bar * F;   % (N x nx) — entrée du bloc

% Hessien symétrisé
H_qp = G' * Q_bar * G + R_bar;
H_qp = (H_qp + H_qp') / 2;

% Contraintes : u_min <= U <= u_max
lb = u_min * ones(N, 1);
ub = u_max * ones(N, 1);