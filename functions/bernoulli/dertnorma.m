function dtnormav2=dertnorma(CIx,CIy,t,deltat0,velocity)
v2=velocity(t+deltat0,CIx,CIy);
v1=velocity(t,CIx,CIy);
dtnormav2(:,:)=(norma(v2)-norma(v1))/(deltat0);
end
