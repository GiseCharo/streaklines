function F=inter_bernoulli(t01,t02,t,CIx,CIy,velo,deltat0,x0,nu,dX)
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

%%computation of the eulerian part
tic
[H1,H2]=prebernoullieulerian(CIx,CIy,t,deltat0,nu,dX,velo);
toc
%%computation of the lagrangian part
t0=t01:deltat0:t02;
stline=zeros(length(t0),2);
tic
for i=1:length(t0)
    stline(i,:)=phi(x0,t0(i),t,velo);
end
toc
[H1streak,H2streak]=interpstreak(stline,H1,H2,CIx,CIy,t0);
Dtiphi=dtiphistreak(t0,x0,deltat0,t,velo);

integH1=(sum(H1streak(:)))*deltat0;
H2t=H2streak.*Dtiphi;
integH2t=(sum(H2t(:)))*deltat0;
F=integH1+integH2t;
end


