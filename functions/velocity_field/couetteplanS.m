function dx=couetteplanS(t,x)
%%couette plan stationary
%%U la moitié de la différence de vitesse entre les deux parois ;
%%h moitié de la distance entre les deux parois

h=1;
U=1;
dx = zeros(2,1); 
dx(1)=(U/h)*x(2);
dx(2)=0;
end