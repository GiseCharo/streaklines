function dx=couetteplanNS(t,x)

%%couette plan
%%caso analitico estacionario couette plan
%%U la moitié de la différence de vitesse entre les deux parois ;
%%h moitié de la distance entre les deux parois
%h=1;
%U0=1;

%U(:,:,1)=(U0/h)*CIy;
%U(:,:,2)=0;
dx = zeros(2,1);  
%%couette plan non stationary
h=1;
U0=1;
nu=1;
omega=2*nu*((2*pi/h)^2);

dx(1)=(U0/sinh(2*pi))*(-sin((2*pi/h)*x(2))*cosh((2*pi/h)*x(2))*sin(omega*t)...
    +cos((2*pi/h)*x(2))*sinh((2*pi/h)*x(2))*cos(omega*t));
dx(2)=0;
 
end
