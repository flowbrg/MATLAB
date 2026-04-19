function F = ppfil(Ba,Aa,Tf,visu);

%
% F = ppfil(Ba,Aa,Tf,fig);
% Ba	 : numérateur du modèle augmenté
% Aa    : dénominateur du modèle augmenté
% Tf   : horizon de filtrage
% F    : Polynome de filtrage
% visu : n° de la figure pour la visualisation des pôles 
%			(absence de visualisation par défaut)
%
% Réalise le placement des pôles de filtrage en continu
% (projection dans le demi plan à gauche de -1/Tf)
% 
% Ph Chevrel Sept. 2000


r=roots(Aa);
Tm=1/max([abs(r);1/Tf]);


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
v1=[-1/Tf -1/Tf];
v2=[-1/Tm 1/Tm];
plot(v1,v2,'y');
hold on
	end

for i=1:length(Aa)-1;
	
%
% 	Affichage des racines de Aa(s)
%
	
		if nargin>3
	plot(real(r(i)),imag(r(i)),'og')
	hold on
		end


	if real	(r(i))>0;
	r(i)=-r(i);
	else r(i)=r(i);
	end
	if real(r(i))>-1/Tf;
	r(i)=-1/Tf+sqrt(-1)*imag(r(i));
	else r(i)=r(i);
	end
%
% 	Affichage des racines de F(s)

		if nargin>3
	plot(real(r(i)),imag(r(i)),'+r')
	hold on
		end
	end
    %end

F=real(poly(diag(r)));
zoom on
if nargin>3
pause
end

