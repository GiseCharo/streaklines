function dx=couetteplanSmatrix(t,CIx,CIy)
%%couette plan stationary over a grid of points
%%U la moitié de la différence de vitesse entre les deux parois ;
%%h moitié de la distance entre les deux parois

h=1;
U=1;

dx(:,:,1)=(U/h)*CIy;
dx(:,:,2)=0;
end