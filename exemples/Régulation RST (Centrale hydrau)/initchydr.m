% initchydr.m
%
% Ce programme iniialise les paramètres du modèle de la centrale hydroélectrique


% Modèle et paramètres de la centrale
% ===================================
Ta=10; Tw=1;
TI=0.3; kw=0.5;
[a,b,c,d]=linmod('chydr_s');
centrale=ss(a,b,c,d);
G=tf(centrale); zpk(G)
[ztbo,pbo,k] = ss2zp(a,b,c,d,1)
usat=100;			% saturation ouverture vanne

