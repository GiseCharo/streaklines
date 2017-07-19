function  phi=phi(p0,t0,t,velocity)
%%computation if the flow given an initial time and a point
options = odeset('RelTol',1e-8,'AbsTol',1e-10);
if t>t0
    [T,sol]=ode45(eval(['@' velocity]), [t0 t],p0,options);
    phi=sol(end,:);
elseif t==t0
    phi=p0;
end

end
