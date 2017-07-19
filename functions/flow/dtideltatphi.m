function dtideltatphi=dtideltatphi(x0,t0,deltat0,t,i,velocity)
%%%calcula la derivada de phi tilde en t0
sol2=phi(x0,t0+deltat0,t0+deltat0+t-(i+1)*deltat0,velocity);
sol1=phi(x0,t0,t0+t-(i+1)*deltat0,velocity);

dtideltatphi=(sol2-sol1)/deltat0;
end
