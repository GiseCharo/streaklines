function dx=vortex(t,CIx,CIy)
%%compute the velocity field of a single vortex

dx(:,:,1) = -(CIy-0.5);
dx(:,:,2) = (CIx-0.5);

end