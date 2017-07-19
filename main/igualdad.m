clear all;close all
%%compute the tangeant of the streakline (Dtophi) and his correction (Dtiphi)
%%also compute the equality Dtiphi=Dtophi+Vlag
%%and computate the norm of Dtiphi
%% x0 is the point of injection, deltat0 the delta of time of the injection
%%t0 is  the time interval of injection and t is the last time of injection
x0(1,1)=1;x0(1,2)=0.5;
deltat0=0.01;t=100; t0=0:deltat0:t-deltat0;

Dtophi=zeros(length(t0)-1,2);
Dtiphi=zeros(length(t0)-1,2);
Vlag=zeros(length(t0)-1,2);

for i=1:length(t0)-1
    Dtophi(i,:)=dtophi(x0,t0(i),t,deltat0,velocity);
    Dtiphi(i,:)=dtiphi(x0,t0(i),deltat0,t-t0(i),velocity);
    Vlag(i,:)=lagrangianvelocity(x0,t0(i),t,velocity);
end


funcion=Dtophi +Vlag;

figure()
plot(Dtiphi(:,1),'*')
hold on
plot(funcion(:,1),'*r')
hleg1 = legend('dtiphi','dtophi +vlag');
title('coordenada x')

figure()
plot(Dtiphi(:,2),'*')
hold on
plot(funcion(:,2),'*r')
hleg2 = legend('dtiphi','dtophi +vlag');
title('coordenada y')

%%%computation of the norm of Dtiphi
nor=Dtiphi(:,1).*Dtiphi(:,1)+Dtiphi(:,2).*Dtiphi(:,2);
figure()
plot(t0(1:end-1),nor,'*')
xlabel('ti')
ylabel('norma 2 de Dtiphi')
title('Dtiphi(ti), t fijo=20, detat0min=0.0001,P=(1,0.5)')
