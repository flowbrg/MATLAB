% gramt.m
% function GT=gramt(a,b,T)
% Programme permettant le calcul du grammien transitoire
% par integration numerique
%
% algorithme deduit de Anderson p348
%
% On resout l'equation de Lyapunov : 
%	d/dt(G) = a*G + G*a' + b*b'
% en operant la vectorisation suivante
% vec(d/dtG)=[kron(I,a) + kron(a,I))*vec(G) + vec(b*b')
%
% exemple : T=1; wn=1; z=0.3; [a,b,c,d]=ord2(wn,z); GT=gramt(a,b,T),

function GT=gramt(a,b,T)

[nl,nc]=size(a);
I1=eye(nl,nc);
A=kron(I1,a)+kron(a,I1);
B=vec(b*b');

t=0:T/100:T;
[Y,X] = step(A,B,eye(size(A)),zeros(nl*nc,1),1,t);
YT=Y(length(Y),:); YT=YT';
GT=unvec(YT,nl);



