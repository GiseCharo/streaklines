function dtiphi=dtiphi(x0,t0,deltat0,deltat,velocity)
% compute the derivative of phi tilde in ti
%

if deltat ~= 0
    sol2=phi(x0,t0+deltat0,t0+deltat0+deltat,velocity);
    sol1=phi(x0,t0,t0+deltat,velocity);
    dtiphi=(sol2-sol1)/deltat0;
else
    dtiphi=0;
end
end
