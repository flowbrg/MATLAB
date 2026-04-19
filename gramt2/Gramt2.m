% gramt.m
% function GT=gramt(a,b,T)
% Programme permettant le calcul du grammien transitoire
% par integration numerique
%
% On resout l'equation de Lyapunov : 
%	d/dt(G) = a*G + G*a' + b*b'
% en operant la vectorisation suivante
% vec(d/dtG)=[kron(I,a) + kron(a,I))*vec(G) + vec(b*b')
%
% exemple : T=1; wn=1; z=0.3; [a,b,c,d]=ord2(wn,z); GT=gramt(a,b,T),
%
% cf. P. Chevrel

function GT=gramt2(A,B,T,algo)

n=length(A);
M_inter=zeros(2*n);
M_inter(n+1:2*n,n+1:2*n)=A';
M_inter(1:n,1:n)=-A;
M_inter(1:n,n+1:2*n)=B*B';

frob_m=norm(M_inter*T,'fro');

J_opt=0;
while frob_m/(2^J_opt)>0.5
   J_opt=J_opt+1;
end

t_opt=T/(2^J_opt);

if algo==0
   F0=expm(t_opt*A');
   m_out=expm(M_inter*t_opt);
elseif algo==1
   F0=expm1(t_opt*A');
   m_out=expm1(M_inter*t_opt);
end

Q0=m_out(1:n,n+1:2*n);
Q0=F0'*Q0;

Fk=F0;
Qk=Q0;

for k=0:J_opt-1
   Qk=Qk+Fk'*Qk*Fk;
   Fk=Fk*Fk;
end

GT=Qk;
