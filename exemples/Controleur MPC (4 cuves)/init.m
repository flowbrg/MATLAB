clc; clear; close all

%% 2.1 Modélisation du système
% Initialisation des paramètres
a1 = 1.31e-4;
a2 = 1.51e-4;
a3 = 9.27e-5;
a4 = 8.82e-5;
ga = 0.3;
gb = 0.4;
Sc = 0.06;
g = 9.81;

% Points de reference
h1ref = 0.8;
h2ref = 0.6;

qbref = ((1-ga)*a1*sqrt(2*g*h1ref)-ga*a2*sqrt(2*g*h2ref))/(1-ga-gb);
qaref = (-gb*a1*sqrt(2*g*h1ref)+(1-gb)*a2*sqrt(2*g*h2ref))/(1-ga-gb);
h3ref = ((1-gb)*qbref)^2/(2*g*a3^2);
h4ref = ((1-ga)*qaref)^2/(2*g*a4^2);

% Jacobiennes analytiques
Alin = (1/Sc)*sqrt(g/2).*[-a1/sqrt(h1ref) 0 a3/sqrt(h3ref) 0;
    0 -a2/sqrt(h2ref) 0 a4/sqrt(h4ref);
    0 0 -a3/sqrt(h3ref) 0;
    0 0 0 -a4/sqrt(h4ref)];

Blin = (1/(3600*Sc)).*[ga 0;
    0 gb;
    0 1-gb;
    1-ga 0];

C = eye(4);
D = zeros(4,2);
rank(ctrb(Alin, Blin)); % =4 donc commandable

%% 2.2 Paramètres du modèle

% 2.2.1 Qa_reel = (1+D)*Qa, on peut modéliser l'imprecision par une
% Création du système en représentation d'état
sys = ss(Alin, Blin, C, D);

% Paramètres de simulation
t_end = 3000;
dt    = 5;
t     = (0:dt:t_end)';
N     = length(t);

% Amplitudes des échelons
u1_amp = qaref*3600;
u2_amp = qbref*3600;

% Bruit blanc ±2%
noise1 = u1_amp * 0.02 * (2*rand(N,1) - 1);
noise2 = u2_amp * 0.02 * (2*rand(N,1) - 1);

% Signaux d'entrée : échelon + bruit
% u = [(u1_amp + noise1), (u2_amp + noise2)];
u = [u1_amp*ones(N,1), u2_amp*ones(N,1)];

% Simulation avec lsim
x0 = 0.5*ones(4,1);
[y, t_out, x] = lsim(sys, u, t, x0);

%  Figure 
state_labels = {'x_1 (h_1)', 'x_2 (h_2)', 'x_3 (h_3)', 'x_4 (h_4)'};
colors       = [0.13 0.47 0.71;   % bleu
                0.90 0.33 0.23;   % rouge
                0.20 0.63 0.17;   % vert
                0.60 0.31 0.64];  % violet

figure('Name','Réponse indicielle – États du système', ...
       'Color','w', 'Position',[100 80 1000 700]);

for i = 1:4
    subplot(2,2,i);
    plot(t_out, x(:,i), 'Color', colors(i,:), 'LineWidth', 1.6);
    xlabel('Temps (s)',  'FontSize', 11);
    ylabel(['$' state_labels{i} '$'], 'Interpreter','latex', 'FontSize', 12);
    title(['État ' num2str(i) ' : ' state_labels{i}], 'FontSize', 12);
    grid on;  box on;
    xlim([0 t_end]);
end

%sgtitle({'Réponse à un échelon – Représentation d''état', ...
%         sprintf('u_1 = %.2f  (±2%% bruit),   u_2 = %.2f  (±2%% bruit)', ...
%                 u1_amp, u2_amp)}, ...
%        'FontSize', 13, 'FontWeight', 'bold');
sgtitle({'Réponse à un échelon – Représentation d''état', ...
         sprintf('u_1 = %.2f,   u_2 = %.2f,   sans bruit blanc', ...
                 u1_amp, u2_amp)}, ...
        'FontSize', 13, 'FontWeight', 'bold');

%% 2.2.2
epsilon    = 0.05;
a1_reel    = a1*(1+epsilon);
qbref_reel = ((1-ga)*a1_reel*sqrt(2*g*h1ref)-ga*a2*sqrt(2*g*h2ref))/(1-ga-gb);
qaref_reel = (-gb*a1_reel*sqrt(2*g*h1ref)+(1-gb)*a2*sqrt(2*g*h2ref))/(1-ga-gb);
h3ref_reel = ((1-gb)*qbref)^2/(2*g*a3^2);
h4ref_reel = ((1-ga)*qaref)^2/(2*g*a4^2);
A_reel = (1/Sc)*sqrt(g/2).*[-a1_reel/sqrt(h1ref) 0 a3/sqrt(h3ref) 0;
    0 -a2/sqrt(h2ref) 0 a4/sqrt(h4ref);
    0 0 -a3/sqrt(h3ref) 0;
    0 0 0 -a4/sqrt(h4ref)];
B_reel = (1/Sc).*[ga 0;
    0 gb;
    0 1-gb;
    1-ga 0];
rank(ctrb(A_reel, B_reel)); % =4 toujours commandable
svd(Alin)
svd(A_reel) % La 4eme valeur passe de 2.7e-3 à 2e-4
% Le système devient mal conditionné
% système devient très sensible aux perturbations / bruit
% le contrôle risque de devenir inefficace et/ou instable numériquement

eig(Alin)
eig(A_reel) % Le premier mode passe de 5.4e-3 à 3e-4

% 2.2.3 On peut faire varier gamma_a et gamma_b pour obtenir de nouvelles
% valeurs d'equilibre, telles que la hauteur du bassin concerné soit en 
% dessous de la fuite.

%% Synthèse LQR

xref = [h1ref; h2ref; h3ref; h4ref]
uref = [qaref; qbref]*3600;

xmax = 1.36*ones(4,1);
umax = [3.26; 4];

Q = diag(1 ./ (xmax.^2)); 
R = diag(1 ./ (umax.^2));

[K,~,~] = lqr(Alin, Blin, Q, R);
