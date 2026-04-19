clc;
clear;
close all;


% QTP Identification
% Paramètres du système et du point de fonctionnement
h1max = 1.36; % [m]
h2max = 1.36; % [m]
h3max = 1.30; % [m]
h4max = 1.30; % [m]
hmin  = 0.20; % [m]

qamax = 3.26/3600; % [m3/s]
qbmax = 4/3600; % [m3/s]
qmin  = 0; % [m3/s]

a1 = 1.31e-4; % [m2]
a2 = 1.51e-4; % [m2]
a3 = 9.27e-5; % [m2]
a4 = 8.82e-5; % [m2]

S  = 0.06; % [m2]

ga = 0.70;
gb = 0.60;

h1o = 0.65; % [m]
h2o = 0.66; % [m]
h3o = 0.65; % [m]
h4o = 0.66; % [m]

qao = 1.63/3600; % [m3/s]
qbo = 2/3600; % [m3/s]

g   = 9.81; % [m/s2]

Ts  = 25; % Sampling time [s]

T1  = S*sqrt(2*h1o/g)/a1; % 1ere constante de temps
T2  = S*sqrt(2*h2o/g)/a2; % 2nde constante de temps
T3  = S*sqrt(2*h3o/g)/a3; % 3eme constante de temps
T4  = S*sqrt(2*h4o/g)/a4; % 4eme constante de temps



% Génération des signaux de commande pour l'identification par deux
% démarches (gensig ou idinput)

Nsamples = 2e3;
% Méthode 1 par gensig
Tf = Nsamples*Ts-1;
[u_qa] = gensig("square", randi([1e2,2e3]),Tf,Ts);
[u_qb] = gensig("square", randi([1e2,2e3]),Tf,Ts);

for k= 1:10
    u_qa = u_qa + randi([0,10],1)*gensig("square",randi([1e2,2e3]),Tf,Ts);
    u_qb = u_qb + randi([0,10],1)*gensig("square",randi([1e2,2e3]),Tf,Ts);
end

u_qa = 7e-4*u_qa/max(u_qa);
u_qb = 0.8e-3*u_qb/max(u_qb);















% Méthode 2 par idinput
u_qa = iddata([],u_qa,Ts);
u_qb = iddata([],u_qb,Ts);


vu_qa.time = u_qa.SamplingInstants;
vu_qa.signals.values = u_qa.InputData;
vu_qb.time = u_qb.SamplingInstants;
vu_qb.signals.values = u_qb.InputData;


% Sorties simulées du modèle non-linéaire
% Open Sim4Tanks
s = sim('Sim4Tanks');
Output = [s.h1.time,s.h1.data,s.h2.time,s.h2.data,];

u = [u_qa.InputData, u_qb.InputData];
y = [Output(:,2),Output(:,4)];

data = iddata(y,u,Ts);
figure, plot(data)
% Data pre-processing
data = detrend(data);
figure, plot(data)





%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%

% Partage des données entre données d'identification et de validation
datai = data(1:Nsamples*3/4);
datav = data(Nsamples*3/4+1 : end);

% Identification par procest (Illustration de l'identification non
% paramétrique) Algorithme d'optimisation nonlinéaire TP fminvin solver
% Introduction d'une connaissance physique

opt = procestOptions('InitialCondition', 'estimate','SearchMethod','lm');
[sysproc, offset,lc] = procest(data,{'P1','P2';'P2','P1'},opt)
figure, compare(datav, sysproc)









% Identification par ARX (Illustration de l'identification paramétrique par
% moindre carrés)

opt = armaxOptions('Focus', 'simulation')

na = 2;
nb = 2;
nk = 1;

arx221 = arx(datai, [na*ones(2,2) nb*ones(2,2) nk*ones(2,2)], opt);
figure, compare(datav, arx221)







% Modèle ARMAX
armax_ = armax(datai, arx221, opt)
figure, compare(datav, arx221, armax_)










% Comparaison avec le modèle linéarisé tangent autour du point de
% fonctionnement
A = [-1/T1 0 1/T3 0;
    0 -1/T2 0 1/T4;
    0 0 -1/T3 0;
    0 0 0 -1/T4];
B = (1/S)*[ga 0;
    0 gb;
    0 (1-gb);
    (1-ga) 0];
C = [1 0 0 0; 0 1 0 0];
D = zeros(2);
syslin = ss(A,B,C,D);
figure, compare(datav,syslin)



% Identification d'un modèle d'état non structuré (ssest)
% Démarche sous jacente méthode des sous espaces et pem
% n4sid(datai, 1:6)
nx = 4;
sysid = ssest(datai,nx,'focus','simulation','DisturbanceModel','none','InitialState','estimate');
x0 = [h1o h2o h3o h4o];
% opt = CompareOptions('InitialCondition',x0)
figure, compare(datav,sysid)



% Identification forme d'état fixée
K = zeros(4,2);
Te = 0;
x0 = [h1o h2o h3o h4o];
init_sys = idss(A,B,C,D,K,x0,Te);

init_sys.Structure.A.Free = (A~=0);
init_sys.Structure.B.Free = (B~=0);
init_sys.Structure.C.Free = false;
init_sys.Structure.D.Free = false;
init_sys.Structure.K.Free = false;

sysid_struct = ssest(datai, init_sys, 'focus', 'simulation','DisturbanceModel','none','InitialState','estimate');
figure, compare(datai,sysid_struct)









function dh = Four_TanksFun(u)
a1 = 1.31e-4;
a2 = 1.51e-4;
a3 = 9.27e-5;
a4 = 8.82e-5;
S = 0.06;
ga = 0.70;
gb = 0.60;
g = 9.81;
h1 = u(1);
h2 = u(2);
h3 = u(3);
h4 = u(4);
qa = u(5);
qb = u(6);
dh1 = (-a1/S)*((2*g*h1)^(1/2)) + (a3/S)*((2*g*h3)^(1/2)) + (ga/S)*qa;
dh2 = (-a2/S)*((2*g*h2)^(1/2)) + (a4/S)*((2*g*h4)^(1/2)) + (gb/S)*qb;
dh3 = (-a3/S)*((2*g*h3)^(1/2)) + ((1-gb)/S)*qb;
dh4 = (-a4/S)*((2*g*h4)^(1/2)) + ((1-ga)/S)*qa;
dh = [dh1; dh2; dh3; dh4];
end