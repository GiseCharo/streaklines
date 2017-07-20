function dtnormav2=dertnorm2(CIx,CIy,t,deltat0,velocity)
v2=velocity(t+deltat0,CIx,CIy);
v1=velocity(t,CIx,CIy);
dtnormav2(:,:)=(norm2(v2)-norm2(v1))/(deltat0);
end
