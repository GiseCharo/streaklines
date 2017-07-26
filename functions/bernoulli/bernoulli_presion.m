function pt02=bernoulli_presion(t01,t02,rho,t,pt01,CIx,CIy,velocity,deltat0,x0,nu,dX)

F=inter_bernoulli(t01,t02,t,CIx,CIy,velocity,deltat0,x0,nu,dX);
vlagt02=lagrangianvelocity(x0,t02,t,velocity);
vlagt01=lagrangianvelocity(x0,t01,t,velocity);
n2=norm2(vlagt02);n1=norm2(vlagt01);

pt02=rho*(F-((n2-n1)/2)+(1/rho)*pt01);
end