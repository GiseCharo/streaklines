function vlag=lagrangianvelocity(x0,t0,t,velocity)
%compute the lagrangian velocity
flux=phi(x0,t0,t,velocity);
vlag=velocity(t,flux(1,1),flux(1,2));
end
