% TP_Systune_head_disk_Matlab
% This example uses a 9th-order model of the head-disk assembly (HDA) in a hard-disk drive. This model captures the first few flexible modes in the HDA.
% On considère ici la problématique de l’asservissement de la position de l'ensemble tête-disque d'un lecteur de disque dur(HDA pour Hard Disk Assembly) .
% Le modèle utilisé, d'ordre 9, capture les premiers modes flexibles du HDA (cf. Mathworks).
% La structure de contrôle utlisée se compose d'un contrôleur PI et d'un filtre passe-bas dans la chaîne de retour. 
% Pour le détail des spécifications, cf. TP ad hoc, Ph. Chevrel, UE PERACT/cours Commande Multiobjectif, TAF ASCy, IMT Atlantique campus Nantes
% L'exercice est ici traité sous Matlab, au travers de l'application systune

load HDA_model G    %load rctExamples G
% load C:\Users\xxxxx\Documents\MATLAB\Examples\systuneWorkflowExample\rctExamples G
figure(1), bode(G), grid, zpk(G), title('modèle HDA')

% You can use the tunablePID object to parameterize the PI block
% To parameterize the lowpass filter F(s), create a tunable real parameter a and construct a first-order transfer function with numerator a and denominator s+a:
C0 = tunablePID('C','pi');  % tunable PI
a = realp('a',1);    % filter coefficient
F0 = tf(a,[1 a]);    % filter parameterized by a
 
% Next build a closed-loop model of the feedback loop. 
% To facilitate open-loop analysis and specify closed-loop requirements such as desired perf or robustness margins, 
% add analysis point : cf. Figure 2 du TP ad hoc, cours Commande Multiobjectif
% Use feedback to build a model of the closed-loop transfer from reference r to head position y
% The result T0 is a generalized state-space model (genss) that depends on the tunable elements C and F

Gu = AnalysisPoint('u');
Gv = AnalysisPoint('v');
Ge = AnalysisPoint('e');
Gz = AnalysisPoint('z');
Gy1 = AnalysisPoint('y1');

T0 = feedback(Gz*G*Gu*Gv*C0*Ge, F0*Gy1);  % closed-loop transfer from r to y (GC/(1+GCF))
T0.InputName = 'r';
T0.OutputName = 'y';
getPoints(T0) % vérification des points d'analyse

% Définition des contraintes
p=tf('s');

gamma_rob1 = 1.2;
gamma_rob2 = 1e3;

% ReqSoft= TuningGoal.Gain('r','e',1e-3*p);               % contrainte soft sur Ter : |Ter(jw)|<|jw|
% ReqHard1= TuningGoal.Gain('z','y',1.2);                 % contrainte hard sur S (sensibilité en sortie)
ReqSoft= TuningGoal.Gain('v','u',1e-3*p);                 % contrainte soft sur Ter : |Ter(jw)|<|jw|
ReqHard1= TuningGoal.Gain('v','u',gamma_rob1);                   % contrainte hard sur S (sensibilité en sortie)
% ReqHard2= TuningGoal.WeightedGain('u', 'v', p/gamma_rob2, 1);   % contrainte hard sur T (sensibilité complementaire; marge dynamique)
ReqHard2= TuningGoal.WeightedGain('z', 'z', p/gamma_rob2, 1);   % contrainte hard sur T (sensibilité complementaire; marge dynamique)
ReqHard3= TuningGoal.WeightedVariance('y1','u',1/0.05,1);  % contrainte hard sur Sub (sensibilité de la commande); pondération gauche (1) et droite (0.05)
ReqHard4= TuningGoal.Gain('r','y',1.35);                  % contrainte hard sur Tyr

ReqHard = [ReqHard1, ReqHard2, ReqHard3, ReqHard4];

tic
options = systuneOptions('MaxIter', 120000,'Display','sub','RandomStart',4);
[T,fSoft,fhard,info] = systune(T0,ReqSoft,ReqHard,options);
duration = toc;

fSoft
fhard

figure(2);
viewSpec(ReqSoft), title('gabarit critere')

figure(3);
stepplot(G), title('reponse indicielle en boucle ouverte : y/u')

figure(4);
stepplot(minreal(T)), title('Reponse indicielle de y/r')

figure(5);
bode(T), title('Bode y/r')

%showTunable(T) 
getBlockValue(T,'C')
%%
Kp=T.Blocks.C.Kp.Value;
Ki=T.Blocks.C.Ki.Value;
a_opt=T.Blocks.a.Value;

% F_opt = tf(a_opt,[1 a_opt]); zpk(F_opt) % Filtre optimal
% C_opt = tf([Kp Ki],[1 0]); zpk(C_opt)  % PI optimal
F_opt = balreal(tf(a_opt,[1 a_opt])); zpk(F_opt) % Filtre optimal
C_opt = balreal(tf([Kp Ki],[1 0])); zpk(C_opt)  % PI optimal
%%
G=balreal(G);
L=minreal(G*C_opt*F_opt);                    % transfert de boucle
S_opt=balreal(minreal(inv(1+L)));
S_opt2=balreal(minreal(feedback(1,L)));      % fonction de sensibilité (en sortie)
S_comp_opt=1-S_opt;                          % fonction de sensibilité complémentaire
S_comp_opt2=feedback(L,1); norm(S_comp_opt-S_comp_opt2)
Tub=minreal(F_opt*C_opt/(1+G*C_opt*F_opt));  % fonction de sensibilité de la commande
% Tub=minreal(F_opt*C_opt/(1+G*C_opt*F_opt));  % fonction de sensibilité de la commande
disp('norm(p*S_comp_opt,inf)'), norm(p*S_comp_opt,inf)

%%
disp('marges de robustesse'), 
allmargin(L),  
[Gm,Pm,Wcg,Wcp]=margin(L), Mr=Pm*pi/180/Wcp, Mm=1/norm(1/(1+L),'inf'), Tdyn=1/(norm(p*minreal(L/(1+L)),'inf'))
figure(6),
bode(S_opt), grid, title('Bode sensibilité y/d')
figure(7)
bode(S_comp_opt), grid, title('Bode sensibilité complementaire')
figure(8),
bode(Tub), grid, title('Bode sensibilité de la commande : u/b')

%%
figure(9);
viewSpec(ReqSoft,T), title('critere et son gabarit')
figure(10);
viewSpec([ReqHard1 ReqHard2],T), % title('contraintes Hinf sur S et (1-S) et gabarits')
disp('ReqHard3 porte sur la norme H2 de la sensibilité au bruit de mesure qui doit être inférieure à 0.05 : '), norm(minreal(Tub),2)
%viewGoal(ReqHard1,T);

% tracés complémentaires
figure(11),bode(minreal(T),minreal(S_comp_opt)), grid, title('Tyr et 1-S')
figure(12), bode(minreal(1-S_opt),(S_opt)), grid, title('1-S et S')
figure(13), nyquist(L), grid,
    