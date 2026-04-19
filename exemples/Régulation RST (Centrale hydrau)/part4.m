% Fonction de transfert
G1 = tf([1],[TI 1]);
G2 = tf([kw -Tw 1], [kw 0.5*Tw 1]);
G3 = tf([1], [Ta 0]);

G = G1*G2*G3
[zeros_G, poles_G, gain_G]=zpkdata(G,'v')


Tf=1;
Tc=0.3;

B = [kw -Tw 1]
A = [Ta*(TI*kw) Ta*(0.5*Tw*TI+kw) Ta*(TI+0.5*Tw) Ta 0]
Aa = [A 0]

C = ppcom(B,A,Tc)
F = ppfil(B,Aa,Tf)
Abf = conv(C,F)
[Sl,R,ASBR] = bezout(B,Aa,Abf);

R

S = [Sl 0]

T = F*C(end)/B(end)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L = tf(R,S)*tf(B,A);

figure(2);
nyquist(L);
hold off;
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