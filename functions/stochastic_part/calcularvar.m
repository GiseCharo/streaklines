clear ;close all; clc
%%compute the variance of the lagrangian velocity
x0(1,1)=1;x0(1,2)=0.5;
t=10;deltat0=0.0005;
t0=0:deltat0:t;
vel='DGyreNS';
for k=1:length(t0)
       Vlag(k,:)=lagrangianvelocity(x0,t0(k),t,vel);
      
end

varVlag=varianza(Vlag);