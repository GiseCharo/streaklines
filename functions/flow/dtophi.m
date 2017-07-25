function dtophi=dtophi(x0,t0,t,deltat0,velocity)
% Compute the tangeante of the streakline on the point ....

if t ~= t0
    sol2=phi(x0,t0+deltat0,t,velocity);
    sol1=phi(x0,t0,t,velocity);
    dtophi=(sol2-sol1)/deltat0;
else
    dtophi= - velocity_on_pt(velocity,t,x0);
end

end
