function U=couetteplanNSmatrix(t,CIx,CIy)

%%couette plan computated over a matrix of points

%%couette plan non stationary
h=1;
U0=1;
nu=1;
omega=2*nu*((2*pi/h)^2);

U(:,:,1)=(U0/sinh(2*pi))*(-sin((2*pi/h)*CIy).*cosh((2*pi/h)*CIy).*sin(omega*t)...
    +cos((2*pi/h)*CIy).*sinh((2*pi/h)*CIy).*cos(omega*t));
U(:,:,2)=0;
 
end