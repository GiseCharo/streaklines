clear; clc;
%%compute the function H (eulerian part of bernoulli equation) for all
%%the points in the domain.
%%H1=(1/2)*dtnormav2-nu*<v,lapv>;
%%H2=vxw-dtv+nu*lapv;
%%Then it computes the lagrangian part, and 
%%then it does the integration along a part of a streakline.
%%[t01 t02] is the integration time interval, t the last time
%%of the injection
%%deltat0=time step of the injection, x0 point of injection
%%CIx,CIy the domain of the problem

t01=0;t02=1;t=2;
dX=1e-4*[1 1];
X=-0.05:dX(1):0.05;
Y=-0.05:dX(2):0.05;
[CIx,CIy]=ndgrid(X,Y);
%%velocity field
vel=@couetteplanNS;
%vel=@couetteplanS;
U0=1;h=1;nu=1;omega=2*nu*((2*pi/h)^2);

deltat0=(pi/omega)/1000;%% dt<<pi/omega

%%computation of the eulerian part
[H1,H2]=prebernoullieulerian(CIx,CIy,t,deltat0,nu,dX,vel);

%%computation of the lagrangian part
x0(1,1)=0.01;x0(1,2)=0.01;
t0=t01:deltat0:t02;
stline=zeros(length(t0),2);

for i=1:length(t0)
stline(i,:)=phi(x0,t0(i),t,vel);
end
[H1streak,H2streak]=interpstreak(stline,H1,H2,X,Y,t0);
Dtiphi=dtiphistreak(t0,x0,deltat0,t,vel);

integH1=(sum(H1streak(:)))*deltat0;
H2t=H2streak.*Dtiphi;
integH2t=(sum(H2t(:)))*deltat0;

