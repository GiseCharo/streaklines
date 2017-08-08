clear;close all
%%computes the streakline
%%x0 is the point of injection,t0 the first time of injection and t the last,
%%deltat0 is the time interval of injection,


% x0(1,1)=1;x0(1,2)=0.5;
% t0=0;deltat0=0.01;t_final=100;
% % vel=@couetteplanS;
% % vel=@DGyreS;
% % vel=@DGyreNS;
% % vel=@couetteplanNS;

% x0(1,1)=6;x0(1,2)=0;
% t0=0;deltat0=0.1;t=10;
% vel=@fct_wake;

% x0(1,1)=1;x0(1,2)=0.5;
% t0=0;deltat0=0.01;t=15;
% % vel=@couetteplanS;
% % vel=@DGyreS;
% % vel=@DGyreNS;
% %vel=@couetteplanNS;
% %vel=@saddle;
% %vel=@fct_wake_megaRAM;

% t0=0;deltat0=0.01;t=20;
t0=0;deltat0=1e-1;t_final=10;
x0(1,1)=6;x0(1,2)=0;
vel=@fct_wake_megaRAM;

vect_time = t0:deltat0:t_final;
N = length(vect_time);
kurt = [0 nan(1,N-1)];
nb = [0 nan(1,N-1)];
nb_by_time = [0 nan(1,N-1)];
nb_by_time2 = [0 nan(1,N-1)];
subsample=50;
figure(1);
taille_police=12;


folder_res = ['images/' func2str(vel) '/kurtosis_vs_time/'];
mkdir(folder_res)

for j=2:subsample:N
    t_local=vect_time(j);
    vect_t0=t0:deltat0:t_local;
    streakline1=streakline(x0,t0,deltat0,t_local,vel);
    Dt0phi = (streakline1(2:end,:)-streakline1(1:end-1,:))/deltat0;
    Dt0phi(end+1,:) = - velocity_on_pt(vel,t_local,x0);
%     Dt0phi = dt0phistreak(vect_t0,x0,t,vel);
    norm_Dt0phi = sum(Dt0phi .^2,2);
%     norm_Dt0phi = sum(dt0phistreak(vect_t0,x0,t_local,vel) .^2,2);
    kurt(j) = fct_kurtosis(norm_Dt0phi);
    nb(j) = fct_nb_peak(norm_Dt0phi);
    nb_by_time(j) = nb(j)/(t_local-t0);
    nb_by_time2(j) = nb(j)/(t_local-t0)^3;
    
    subplot(3,2,1)
    plot(streakline1(:,1),streakline1(:,2),'o-')
    title('Streakline')
    subplot(3,2,2)
    plot(vect_t0,norm_Dt0phi,'-')
    title('Tangeante norm')
    subplot(3,2,3)
    plot((t0+deltat0):(subsample*deltat0):t_local,...
        kurt(2:subsample:j)-3,'-')
    title('Excess kurtosis')
    subplot(3,2,4)
    plot((t0+deltat0):(subsample*deltat0):t_local,...
        nb(2:subsample:j),'-')
    title('Nb jumps')
    subplot(3,2,5)
    plot((t0+deltat0):(subsample*deltat0):t_local,...
        nb_by_time(2:subsample:j),'-')
    title('Nb jumps by unit of time $t_0$ (s$^{-1}$)',...
        'interpreter','latex',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontSize',taille_police,...
        'FontName','Times')
    subplot(3,2,6)
    plot((t0+deltat0):(subsample*deltat0):t_local,...
        nb_by_time2(2:subsample:j),'-')
    title('Nb jumps by unit of $t_0$ and by unit of $t^2$ (s$^{-2}$)',...
        'interpreter','latex',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontSize',taille_police,...
        'FontName','Times')
    drawnow;    
    
    
    eval( ['print -depsc ' folder_res 't=_' num2str(t_local) '.eps']);
end