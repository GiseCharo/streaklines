function dx=couetteplanNS(t,x)
%%couette plan non stationary
%%the parameters must be changed in the main script, also here
dx = zeros(2,1);  
h=1;U0=1;nu=1;omega=2*nu*((2*pi/h)^2);

dx(1)=(U0/sinh(2*pi))*(-sin((2*pi/h)*x(2))*cosh((2*pi/h)*x(2))*sin(omega*t)...
    +cos((2*pi/h)*x(2))*sinh((2*pi/h)*x(2))*cos(omega*t));
dx(2)=0;
 
end
