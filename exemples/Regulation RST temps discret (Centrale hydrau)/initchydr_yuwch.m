clc; clear; close all;
% initchydr_yuwch.m
%
% Ce programme iniialise les paramètres du modèle de la centrale hydroélectrique


% Modèle et paramètres de la centrale
% ===================================
Ta=10; 
Tw=1;
TI=0.3; 
kw=0.5;

usat=100;			% saturation ouverture vanne

[a,b,c,d]=linmod('chydr_s_yuwch');
b_u = b(:,1);
b_p = b(:,2);
d_u = d(:,1);
d_p = d(:,2);
Gssyuwch=ss(a,b,c,d);

Gyu=tf(Gssyuwch(:,1)); 
G = Gyu

Gywch=tf(Gssyuwch(:,2)); 

load("datahydr_pprTf5.mat");
%% Calcul d'un régulateur RST discret à partir du RST continu
Te = 0.4;

RS      = tf(R,S);             % R/S
RSd     = c2d(RS,Te,'tustin'); % Rd/Sd
[Rd,Sd] = tfdata(RSd,'v');

TS      = tf(T,S);             % T/S
TSd     = c2d(TS,Te,'tustin'); % Td/Sd
[Td,Sd] = tfdata(TSd,'v');

% Polynome de filtrage
[Fd,~]  = tfdata(zpk([0 0 0 0 0],[],1,Te),'v');


% Diagrammes de Bode

figure
bode(RS, RSd, TS, TSd)
legend('R/S continu', 'R/S discret', ...
       'T/S continu', 'T/S discret')


%% Reconception d'un régulateur RST
Te = 0.4;

Cd = poly(exp(Te*roots(C))) % Polynome de commande discret
Fd = poly(exp(Te*roots(F))) % Polynome de filtrage discret

% Action intégrale donc on a S_d(1) = 0 càd S_d(z) = (z-1) * S_L_d(z)
% L'equation de Bezout donne Abf_d(z) = A_d(z) * S_d(z) + ....
% Soit A_d(z)*(z-1)*S_L_d(z); on pose alors A_augmenté_d(z) = (z-1)*A_d(z)

[Bd,Ad] = tfdata(c2d(G,Te,'zoh'),'v');
Aad = conv(Ad,[1 -1]); % A_augmenté_d = (z-1)*A_d

Abf = conv(Cd,Fd);
[SLd,Rd] = bezout(Bd,Aad,Abf);

Sd = conv(SLd,[1 -1]);

% Condition Syyc(1) = 1 donc (B_d(1) * T_d(1))/(C_d(1) * F_d(1)) = 1
% D'où T_d(z) = B_d(1)/C_d(1) * F_d(z)

Td = Fd*sum(Cd)/sum(Bd);

%% Conception d'un régulateur par retour d'état observateur visant des performances équivalentes
Te = 0.4;

% p = tf("s");
% C*inv(p*eye(4)-A)*B(:,1)

% Modele discret
[ad,bd,cd,dd] = c2dm(a,b,c,d,Te,'zoh');
% Gssyuwch_d    = ss(ad,bd,cd,dd,Te);

% Retour d'état
bd_p = bd(:,2); % Matrice B du signal de commande
bd_u = bd(:,1); % Matrice B de la perturbation

Kd   = place(ad,bd_u, roots(Cd)) % K = place(ad,bd(:,1), exp(Te*roots(C)))
% eig(ad-bd_u*K)
% exp(Te*roots(C))
% Valeurs propres identiques aux exp(Te*roots(C))

% Coeffecients correcteurs
Kdc  = inv(cd*inv(eye(4)-ad+bd_u*Kd)*bd_u) 
Kdp  = -inv(cd*inv(eye(4)-ad+bd_u*Kd)*bd_u)*cd*inv(eye(4)-ad+bd_u*Kd)*bd_p

% Observateur

% Modele augmenté discret
aad  = [ad bd_p; zeros(1,4) 1];
bad  = [bd_u; 0];
cad  = [cd 0];
dad  = 0;

Ld   = place(aad', cad', exp(Te*roots(F)))' % Gain de l'observateur
Aobs = aad-Ld*cad;
Bobs = [bad Ld];
Cobs = eye(5);
Dobs = zeros(5,2);