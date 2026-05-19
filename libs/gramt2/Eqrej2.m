% EQREJ2.M  
%
% function [Ta, Ga]=eqrej(A,B,C,n1);
%
% programme permettant la résolution des équations
% du rejet asymptotique :
%
%	A11*Ta - Ta*A22 + B1*Ga = A12
%	C1*Ta 			= C2
%
% paramètres d'entrée :
%	# A=[A11 A12 ; 0 A22]
%	# B=[B1 ; 0]
%	# C=[C1 , C2]
%	# n1 : nombre d'états commandables (taille de A11)
%
% Auteur : P. Chevrel
% date de mise à jour : 10/95
%

function [Ta, Ga]=eqrej2(A,B,C,D,n1);

tol=1e-8;
[n,m]=size(B);
[p,n]=size(C);
n=length(A);
n2=n-n1;
nu=m;

A11=A(1:n1,1:n1);	
A12=A(1:n1,n1+1:n);
A21=A(n1+1:n,1:n1);
A22=A(n1+1:n,n1+1:n);
B1=B(1:n1,:);
B2=B(n1+1:n,:);
C1=C(:,1:n1);
C2=C(:,n1+1:n);


if (max(max(abs(A21))) > tol)|(max(max(abs(B2))) > tol)
  error('le système doit être donné sous forme canonique de commandabilité')
end




[U1,T1]=schur(A22);
[U,A22n]=rsf2csf(U1,T1);

A12n=A12*U;
C2n=C2*U;

for ii=1:1:n2
   S=zeros(n1,1);
   for kk=1:1:ii-1;
      S=S+Ta1(1:n1,kk)*A22n(kk,ii);
   end
   P=[A11-A22n(ii,ii)*eye(n1) B1; C1 D];
   R=[A12n(:,ii)+S;C2n(:,ii)];
   F=inv(P)*R;
   Ta1(1:n1,ii)=F(1:n1,1);
   Ga1(1:nu,ii)=F(n1+1:n1+nu,1);
end
Ta=Ta1(1:n1,1:n2)*U';
Ta=real(Ta);
Ga=Ga1(1:nu,1:n2)*U;
Ga=real(Ga);
   
      