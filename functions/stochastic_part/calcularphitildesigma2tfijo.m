clear; close all;clc
%% Compute the sigma2(deltat0)of the correction of the tangent of the flux
%%when  t  is fix,
x0(1,1)=1;x0(1,2)=0.5;
t=20; deltat0min=0.0001;
nmax=20;
t0=deltat0min:deltat0min:t;
f=zeros(length(t0),2);
vel='DGyreNS';
for i=1:length(t0)
    f(i,:)=phitilde(x0,t0(i),t-t0(i),vel);
end
clear t0;
deltat0=deltat0min:deltat0min:t;
%deltat0=p*deltat0min;
ind=0;
for j=1:length(deltat0)
    p=round(deltat0(j)/deltat0min);
    t0=0:p*deltat0min:t;
    if length(t0)-1>nmax
        ind=ind+1;
        g=zeros(length(t0)-1,2);
        for n=1:length(t0)-2
            g(n,:)=f(p*n+1,:)-f(p*(n-1)+1,:);
        end
        sigma2dtiphitfijo(ind)=varianzabb(g,deltat0(j));
        coorddeltat01(ind)=deltat0(j);
         clear g t0;
    end
end

figure()
plot(coorddeltat01,sigma2dtiphitfijo,'o-')
title('sigma2dtiphi(deltat0), t fijo=20, detat0min=0.005,P=(0.5,0.5)')
xlabel('deltat0')
ylabel('sigma2dtiphi')

