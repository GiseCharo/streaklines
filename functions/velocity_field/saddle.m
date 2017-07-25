function dx=saddle(t,CIx,CIy)
%%compute the velocity field of the velocity field with a saddle point in 
%%(0,0)
omega=pi/4;
lambda=2;
A=0.1;
forcing=A*cos(omega*t);
dx(:,:,1) = lambda*CIx+forcing;
dx(:,:,2) = -lambda*CIy;

end