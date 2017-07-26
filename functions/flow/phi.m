function  phi=phi(p0,t0,t,velocity)
%%computation if the flow given an initial time and a point
%

velocity=@(t,CIx) velocity_on_pt(velocity,t,CIx);
options = odeset('RelTol',1e-8,'AbsTol',1e-10);
% %if t ~= t0
% if t>t0
%     [~,sol]=ode45(velocity, [t0 t],p0,options);
%     phi=sol(end,:);
% elseif t==t0
%     phi=p0;
if t>t0+10*eps
    [~,sol]=ode45(velocity, [t0 t],p0,options);
    phi=sol(end,:);
elseif abs(t-t0)<10*eps
    phi=p0;
else
    error('Backward advection not encoded yet');
end
end
