%% This script advect and plot a grid of points
% It also computes several mixing diagnoses (Mezic, FTLE, Okubo-Weiss,...)
% For the Lagrangian mixing diagnoses (Mezic, FTLE), the initial time t_0
% is fixed but the advection time, tau, increases.
% The boolean parameter back choose between a forward or a Backard
% advection.

init

% Advection backward or forward
%back = true;
back = false;

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
%
% % Non-stationnary double gyre
% t0=0;t=50;
% N= 50;
% vel=@DGyreNS;
% x=linspace(0,2,2*N);
% y=linspace(0,1,N);
% T_plot = 10;
% n_sub = 1;

% % Wake flow
% t0=0;t=20;
% N= 120;
% vel=@fct_wake;
% x=linspace(0,20,2*N);
% y=linspace(-6,6,N);
% T_plot = 100;
% n_sub = 3;

% Wake flow MEGA RAM
% t0=1;t=15;
% vel=@fct_wake_megaRAM;
t0=1;t=35;
vel=@fct_wake_megaRAM_2blocks;
% % N= 120;
% % x=linspace(-20,20,4*N);
% % %x=linspace(0,20,2*N);
% % y=linspace(-6,6,N);
% N = 100;
% x=linspace(5,20,10*N);
% y=linspace(-0.5,0.5,N/2);
% %y=linspace(-1,1,N);
%N = 5;
 N = 20;
% N = 100;
x=linspace(-10,20,7.5*N);
%x=linspace(5,20,7.5*N);
y=linspace(-3,3,3*N);
%y=linspace(-1,1,N);
T_plot = 500;
% T_plot = 100;
n_sub = 1;

% Grid
[X,Y]=ndgrid(x,y);

% Lagrangian advection
if back
    [X,Y] = phi_grid_back(X,Y,t0,t,vel,true,n_sub,T_plot );
else
    [X,Y] = phi_grid(X,Y,t0,t,vel,true,n_sub,T_plot );
end