clear;clc;
%%Compute the covariance of the lagrangian velocity and of the correction of
%%the tangent of the flux
x0(1,1)=1;x0(1,2)=0.5;
t=10;deltat0=0.01;
t0=0:deltat0:t;
Vlag1=zeros(length(t0),2);
Dtiphi1=zeros(length(t0),2);
vel=@DGyreNS;
%vel=@DGyreS;
%vel=@couetteplanNS;
%vel=@couetteplanS;
for i=1:length(t0)
    Vlag1(i,:)=lagrangianvelocity(x0,t0(i),t,vel);
    Dtiphi1(i,:)=dtiphi(x0,t0(i),deltat0,t-t0(i),vel);
end
tam=floor(length(t0)/2);
tau0=t0(1:tam);
covV=covarianza(Vlag1);
covdtiphi=covarianza(Dtiphi1);


figure()
plot(tau0,covV)
title('covV(tau)')
xlabel('tau0')
ylabel('covV')

figure()
plot(tau0,covdtiphi)
title('covdtiphi(tau)')
xlabel('tau0')
ylabel('covdtiphi')
