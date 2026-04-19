function [S,R] = bezout(B,A,C)
%Bezout 
%
%   Résout l'équation de Bezout :
%
%		A(p)S(p)+B(p)R(p)=C(p)
%
%   Le couple de solutions R et S obtenues à partir de A, B et C peut être obtenu par résolution du système
%	 d'équations linéaires :
%
%		MeX=Y
%
% 	 où :
%			* Me est la matrice de Sylvester déduite des coefficients de A et B selon : (n=na)
%			        | a0             b0             |
%			        | a1             b1             |
%			        | .              .              |
%			        | .              .              |
%			        | .              .              |
%			    Me = | an             bn             |
%			        |     an             bn         |
%			        |         an             bn     |
%			        |            an              bn |
%   		* Y est le vecteur déduit des coefficients de Abf
%				 Y = [ c0 c1 ...	ci ...	cnc ]'
%
%	 Bibliographie :    G. C. Goodwin & K. S. Sin., "Adaptive Filtering, Prediction, and Control",
%
%   Philippe Chevrel
%   nov. 2015


nA = length (A) - 1; 
nB = length (B) - 1;
if nB>nA disp('système impropre'), return, end
nC=length(C)-1;
nS=nC-nA;
nR=nC-nS-1;

A; B; C;
Me = zeros (nC+1,nS+nR+2);
for i = 1:nS+1
  Me(i:i+nA,i) = fliplr(A)';				% A=[anA ... an0] => fliplr(A)=[an0 ... anA]
end
for i = 1:nR+1
  Me(i:i+nB,i+nS+1) = fliplr(B)';
end
Me

if rcond (Me) > eps
  SR = Me \ fliplr(C)';
  S = fliplr(SR(1:nS+1)');
  R = fliplr(SR(nS+2:nS+nR+2)');
  Rl=zeros(size(C)); Rl(nC-nR+1:nC+1)=R, 
  AS=conv(A,S);
  BR=conv(B,R); BRl=zeros(size(C));BRl(nC-nR-nB+1:nC+1)=BR;  %BRl=[0 0 BR]; mise a la bonne longueur
  ASBR=AS+BRl;
  disp(['erreur maximum sur la localisation des pôles en boucle fermée  : -->  ',num2str(max(abs(roots(ASBR)-roots(C))))]);
else
  disp('The eliminant matrix is singular.');
  S = [];
  R = [];
end






