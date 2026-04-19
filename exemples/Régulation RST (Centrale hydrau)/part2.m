% Fonction de transfert
G1 = tf([1],[TI 1]);
G2 = tf([kw -Tw 1], [kw 0.5*Tw 1]);
G3 = tf([1], [Ta 0]);

G = G1*G2*G3
[zeros_G, poles_G, gain_G]=zpkdata(G,'v');

disp('========================================')
% Display the transfer functions
disp('Transfer Function G:');
gain_G
zeros_G
poles_G
disp('========================================')

figure(1);
step(G)
figure(2);
bode(G)

% Regulateur
[zeros_K, poles_K, gain_K] = zpkdata(K,'v');

disp('========================================')
% Display the transfer functions
disp('Transfer Function K:');
gain_K
zeros_K
poles_K
disp('========================================')


L = G*K;


% Fonction de sensibilité
Syd = 1/(1+L);
[zeros_Syd, poles_Syd, gain_Syd]=zpkdata(Syd,'v');
disp('========================================')
% Display the transfer functions
disp('Transfer Function K:');
gain_Syd
zeros_Syd
poles_Syd
disp('========================================')

figure(3);
step(Syd);

Sc = 1-Syd;
Sub = K/(1+L);

figure(4);
bode(Syd);
hold on;
bode(Sc);
bode(Sub);
hold off;
legend('Syd(jω) - Sensibilité à la perturbation', 'Sc(jω) - Complémentaire', 'Sub(jω) - Sensibilité de la commande au bruit', 'Location', 'best');
title('');
grid on;

figure(5);
nyquist(L);
title('Diagramme de Nyquist de la chaîne directe')

figure(6);
bode(L);
title('Diagramme de Bode de la chaîne directe')

% Calcul des marges
[Gm, Pm, Wcg, Wcp] = margin(L);
disp('========================================');
disp('Analyse de la châine directe');
disp('========================================');
disp(['Marge de gain      : ', num2str(20*log10(Gm)), ' dB']);
disp(['Marge de phase     : ', num2str(Pm), ' degrés']);

[re, im] = nyquist(L);
distances = sqrt((re(:) + 1).^2 + im(:).^2);
d_min = min(distances);
disp(['Distance min au point critique : ', num2str(d_min)]);
disp(' ');

[mag, ~] = bode(1+L);
Mm = min(mag(:));
disp(['Marge de module (Mm) : ', num2str(20*log10(Mm)), ' dB']);
disp('========================================');
