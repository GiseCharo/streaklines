%% This script advect and plot a grid of points

init

% % Stationnary Couette-plan
% t0=0;t=50;
% N= 50;
% vel=@couetteplanS;
% x=linspace(0,2,2*N);
% y=linspace(0,1,N);
% T_plot = 30;
% n_sub = 3;

% % Non-stationnary Couette-plan
% t0=0;t=50;
% N= 50;
% vel=@couetteplanNS;
% x=linspace(0,2,2*N);
% y=linspace(0,1,N);
% T_plot = 30;
% n_sub = 3;

% % Stationnary double gyre
% t0=0;t=50;
% N= 50;
% vel=@DGyreS;
% x=linspace(0,2,2*N);
% y=linspace(0,1,N);
% T_plot = 10;
% n_sub = 1;

% % Non-stationnary double gyre
% t0=0;t=50;
% N= 50;
% vel=@DGyreNS;
% x=linspace(0,2,2*N);
% y=linspace(0,1,N);
% T_plot = 10;
% n_sub = 1;

% % Wake flow
% t0=0;t=10;
% N= 120;
% vel=@fct_wake;
% x=linspace(0,20,2*N);
% y=linspace(-6,6,N);
% T_plot = 100;
% n_sub = 3;

% Wake flow MEGA RAM
t0=0;t=10;
N= 120;
vel=@fct_wake_megaRAM;
x=linspace(0,20,2*N);
y=linspace(-6,6,N);
T_plot = 100;
n_sub = 3;

% Grid
[X,Y]=ndgrid(x,y);

% Lagrangian advection
[X,Y] = phi_grid(X,Y,t0,t,vel,true,n_sub,T_plot );
