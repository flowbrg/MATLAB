%% INITIALISATION - Moteur synchrone biphase à champ statorique induit
% Cette modelisation est realisee dans le referenciel du stator
% A executer avant de lancer le modele


clear; clc;

%% Parametres geometriques et magnetiques
mu0     = 4*pi*1e-7;    % Permeabilite du vide          [H/m]
R       = 0.1;          % Rayon moyen entrefer          [m]
e       = 1e-2;         % Epaisseur entrefer            [m]
L_ax    = 0.2;          % Longueur axiale               [m]
V_ent   = 2*pi*R*e*L_ax; % Volume entrefer              [m^3]

N_spire = 80;           % Nombre de spires
kb      = 0.9;          % Facteur de bobinage 
Ks      = mu0*N_spire*kb/(2*e); % Facteur bobinage statorique [T/A]
Kr      = Ks;           % Facteur bobinage rotorique    [T/A]
M0      = (V_ent/(2*mu0))*Ks*Kr; % Mutuelle inductance max [H]

%% Parametres electriques rotoriques
Rr      = 2.0;          % Resistance rotorique          [Ohm]
Lr      = 50e-3;        % Inductance propre rotorique   [H]

%% Parametres electriques statoriques
Rs      = 1.5;          % Resistance statorique         [Ohm]
Ls      = 20e-3;        % Inductance propre statorique  [H] 

param_elec = [Rs ; Ls ; Rr ; Lr ; M0];

%% Parametres mecaniques
f_frot  = 1e-3;         % Frottements fluides           [N.m.s/rad]
Gamma_c = 1.0;          % Couple résistant constant     [N.m]
rho     = 7.85e3;       % Masse volumique acier 42CrMo4 [kg/m^3]
m       = pi*R^2*L_ax*rho;  % Masse rotor               [kg]
J       = m*R^2/2;      % Moment d'inertie rotor        [kg.m^2]

%% Parametres source /commande
ws_nom  = 2*pi*5;      % Pulsation nominale 50 Hz      [rad/s]
Usm     = 230;          % Amplitude tension statorique  [V]
Ir_ref  = 4;            % Courant excitation rotorique  [A]
Ur      = Rr * Ir_ref;  % Tension rotorique             [V]

%% Reglage FOC
tau_s   = Ls / Rs;      % Constante de temps electrique [s]
tau_c   = tau_s / 10;   % Constante de temps boucle interne [s]

% PI internes (axes d et q identiques)
Kp_dq   = Ls / tau_c;   % Gain proportionnel            [V/A]
Ki_dq   = Rs / tau_c;   % Gain integral                 [V/(A.s)]
Iq_max = 2 * Gamma_c / (M0 * Ir_ref); % Courant de la phase q saturée [A]

% Boucle externe vitesse
tau_w   = 8 * 10 * tau_c;   % Constante de temps boucle omega [s] hyp : boucle interne 10 fois plus rapide, nombre magique 8
Kp_w    = J / (tau_w * M0 * Ir_ref);  % Gain proportionnel [A.s/rad]
Ki_w    = Kp_w / tau_w; % Gain integral                 [A/rad]

% Reference
id_ref  = 0;            % Courant axe d (pas de magnetisation inutile)
w_ref = ws_nom;     % Consigne vitesse = vitesse synchrone [rad/s]

%% Scenario demarrage
% 'reseau'  : ws = cste = ws_nom (accrochage brutal)
% 'variateur': rampe de 0 a ws_nom
scenario = 'variateur';
t_rampe = 10.0;      % Duree de la rampe                 [s]
kr      = ws_nom / t_rampe; % Pente rampe               [rad/s^2]

%% Conditions initiales
w0      = 0;            % Vitesse initiale              [rad/s]
theta0  = -pi/4;        % Position initiale             [rad/s]
i1_0    = 0;            % Courant initial phase statorique 1 [A]
i2_0    = 0;            % Courant initial phase statorique 2 [A]
Ir_0    = Ur / Rr;      % Courant rotorique initial     [A]
% (rotor pre-excite en courant continu avant t=0)
x_elec_0 = [i1_0; i2_0; Ir_0]; % CI integrateur vectoriel bloc electrique

%% Verification condition d accrochage
Bsm     = Ks * (Usm / sqrt(Rs^2+(ws_nom*Ls)^2)); % Amplitude champ statorique [T]
Brm     = Kr * Ir_ref;                % Amplitude champ rotorique [T]
Gamma_0 = (V_ent / (2*mu0)) * Bsm * Brm; % Couple max   [N.m]
Gamma_charge = Gamma_c + f_frot * ws_nom;    % Couple resistant total [N.m]

fprintf('=== Verification de la condition d accrochage ===\n');
fprintf('Couple max disponible  Gamma_0      = %.2f N.m\n', Gamma_0);
fprintf('Couple resistant total Gamma_charge = %.2f N.m\n', Gamma_charge);
if Gamma_0 > Gamma_charge
    fprintf('OK\n');
else
    fprintf('NOK\n');
end
fprintf('\n');
fprintf('=== Conditions initiales ===\n');
fprintf('ws_nom = %.2f rad/s\n', ws_nom);
fprintf('w0     = %.2f rad/s\n', w0);
fprintf('theta0 = %.2f rad\n', theta0);
fprintf('Ir_ref = %.2f A\n', Ir_ref);
fprintf('Usm    = %.2f V\n', Usm);
fprintf('Ur     = %.2f V\n', Ur);
fprintf('kr     = %.2f rad/s\n', kr);

%% Duree simulation
t_sim   = 30.0;         % Duree totale de simulation    [s]