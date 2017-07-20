function dtv=dertvel(CIx,CIy,t,deltat0,velocity)
v2=velocity(t+deltat0,CIx,CIy);
v1=velocity(t,CIx,CIy);
dtv=(v2-v1)/(deltat0);
end