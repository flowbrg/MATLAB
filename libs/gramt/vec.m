%vec.m
% function V=vec(M)
% Transforme une matrice en vecteur
% 
% M=[m1 m2 m3 ...] avec m1, m2, m3 ... nc vecteurs a nl lignes
% V=[m1 ; m2 ; m3 ; ...] est un vecteur à nl*nc lignes

function V=vec(M)

[nl,nc]=size(M);
V=M(1:nl*nc)';

%V=zeros(nl*nc,1);
%for i=0:nc-1 
%	V(i*nl+1:(i+1)*nl,1)=M(:,i+1); 
%end


