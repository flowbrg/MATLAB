% unvec.m
% function M=unvec(V,nl)
% Transforme un  vecteur colonne en matrice a nl lignes
% 
% V=[m1 ; m2 ; m3 ...] avec m1, m2, m3 ... nc vecteurs a nl lignes
% M=[m1 m2 m3 ...] est une matrice a nl lignes

function M=unvec(V,nl)

[nlV,ncV]=size(V);
if min(nlV,ncV) > 1, error('V est une matrice'), end
if ncV>nlV, V=V'; nlV=ncV; end

nc=round(length(V)/nl);
if abs((length(V)/nl)-nc)>1e-3
	error('la taille du vecteur n"est pas un multiple de nl')
end

for i=1:nc
	M(1:nl,i)=V((i-1)*nl+1:i*nl);
end



