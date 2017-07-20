function dx=DGyreS(t,CIx,CIy)
%%compute the velocity field of the stationary double gyre over a grid
%%of points
epsilon=0.1;
amplitude=0.1;

a = epsilon;
b = 1 - 2*epsilon;
forcing = a*CIx.^2 + b*CIx;

dx(:,:,1) = -pi*amplitude*sin(pi*forcing).*cos(pi*CIy);
dx(:,:,2) = pi*amplitude*cos(pi*forcing).*sin(pi*CIy).*(2*a*CIx + b);

end