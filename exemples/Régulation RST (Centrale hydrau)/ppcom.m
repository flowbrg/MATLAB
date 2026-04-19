function [C] = ppcom(B,A,Tc,visu);

%
% [C] = ppcom(B,A,Tc,fig);
% B    : numérateur du modèle du processus
% A    : dénominateur du modèle du processus
% Tc   : horizon de commande
% C    : Polynome de commande
%
% Réalise le placement des pôles de commande en continu
% selon la méthode PPR [Thèse S. Puren]
% 
% Ph Chevrel Sept. 2000

nA=length(A);
nB=length(poly(roots(B)));     %Calcule length(B), mais en éliminant les "zeros" inutile du pol. B (ex [0 0 1] -> 1 et pas 3
nA_B=nA-nB;
r=-1/Tc*ones(nA-1,1);	
r1=roots(B);
Tm=1/max([abs(r1);1/Tc]);


%
% Tracé des axes 
%

if nargin>3
    figure(visu)
    axis('equal');
    grid on
    zoom on
    hold on
    v1=[-1.5/Tm 0.5/Tm];
    v2=[0 		0];
    plot(v1,v2,'g');
    v1=[0 		0];
    v2=[-1/Tm 	1/Tm];
    plot(v1,v2,'g');
    v1=[-1/Tc -1/Tc];
    v2=[-1/Tm 1/Tm];
    plot(v1,v2,'y');
    hold on
end

	
%
% 	Affichage des racines de B(s)
%
if nargin>3
    plot(real(r1),imag(r1),'og')
    hold on
end
      

% Calcul des racines de C
for i=1:nB-1;
	if real	(r1(i))>0;
	r1(i)=-r1(i);
	else r1(i)=r1(i);
	end
	if real(r1(i))<-1/Tc;
	r1(i)=-1/Tc+sqrt(-1)*imag(r1(i));
	else r1(i)=r1(i);
   end
   r(1:nB-1)=r1;
end   
   
% 	Affichage des racines de C(s)
if nargin>3
	plot(real(r),imag(r),'+r')
	hold on
end

% Calcul de C
C=real(poly(diag(r)));

if nargin>3, zoom on, pause, end

