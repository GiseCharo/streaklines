clear;close all
%%computes the streakline
%%x0 is the point of injection,t0 the first time of injection and t the last,
%%deltat0 is the time interval of injection,


% x0(1,1)=1;x0(1,2)=0.5;
%t0=0;deltat0=0.01;t=5;
% vel=@couetteplanS;
% vel=@DGyreS;
% vel=@DGyreNS;
% vel=@couetteplanNS;

x0(1,1)=6;x0(1,2)=0;
t0=0;deltat0=0.1;t=10;
vel=@fct_wake;

x0(1,1)=1;x0(1,2)=0.5;
t0=0;deltat0=0.01;t=15;
% vel=@couetteplanS;
% vel=@DGyreS;
% vel=@DGyreNS;
%vel=@couetteplanNS;
vel=@saddle;
%vel=@fct_wake_megaRAM;


streakline1=streakline(x0,t0,deltat0,t,vel);

figure()
plot(streakline1(:,1),streakline1(:,2),'o-')
%axis xy;axis equal
%title('streakline, P=(1,0.5), deltat0=0.01, t=30')
%axis([0 2 0 1])
