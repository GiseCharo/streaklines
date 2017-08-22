clear;close all
%% Compute the Q cirterion along a streaklines
init

%% Flow
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

% % % Stationnary double gyre
% x0(1,1)=0.5;x0(1,2)=0.45;
% % x0(1,1)=1;x0(1,2)=1e-3;
% % x0(1,1)=1;x0(1,2)=0.5;
% t0=0.;
% deltat0=1.e-1;
% t=5.;
% deltat=1.e-1;
% % N= 50;
% vel=@DGyreS;
% % x=linspace(0,2,2*N);
% % y=linspace(0,1,N);
% % T_plot = 10;
% % n_sub = 1;
% ax=[0 2 0 1];
% loc_colorbar = 'southoutside';

% Non-stationnary double gyre
x0(1,1)=1;x0(1,2)=0.5;
t0=0.;
deltat0=1.e-2;
%big_deltat0=20;
big_deltat0=10 % GOOD
big_deltat0=5;
t=100.;
deltat=1.e-1;
N= 100.;
vel=@DGyreNS;
x=linspace(0,2,2*N);
y=linspace(0,1,N);
T_plot = 10;
n_sub = 1;
ax=[0 2 0 1];
loc_colorbar = 'southoutside';

% % Single vortex
% x0(1,1)=0.25;x0(1,2)=0.25;
% vel=@vortex;
% t0=0.;
% deltat0=1.e-1;
% t=5.;
% deltat=1.e-1;
% N= 100.;
% x=linspace(0,1,N);
% y=linspace(0,1,N);
% T_plot = 10;
% n_sub = 1;
% ax=[0 1 0 1];
% loc_colorbar = 'southoutside';

% % Wake flow
% t0=0;t=20;
% N= 120;
% vel=@fct_wake;
% x=linspace(0,20,2*N);
% y=linspace(-6,6,N);
% T_plot = 100;
% n_sub = 3;

% % Wake flow MEGA RAM
% t0=0;t=10;
% N= 120;
% vel=@fct_wake_megaRAM;
% x=linspace(0,20,2*N);
% y=linspace(-6,6,N);
% T_plot = 100;
% n_sub = 3;

%% Grid
[X,Y]=ndgrid(x,y);
dX = [X(2,1)-X(1,1) Y(1,2)-Y(1,1)];
grid.dX=dX;
grid.MX=size(X);
grid.x_ref = x;
grid.y_ref = y;
%vect_t0 = t0:deltat0:t;
% vect_t0 = (t0+deltat0):deltat0:(t-deltat0);
vect_t0 = t0:deltat0:(t-big_deltat0);
vect_t0 = vect_t0(end:-1:1);

%% Precomputation

% Streakline
%streakline1=streakline(x0,t0,deltat0,t,vel);
streakline1=streakline_2(x0,vect_t0,t,vel);

% Eulerian velocity
v = vel(t,X,Y);

% Norm of Eulerian velocity
norm_v = sqrt(sum(v.^2,3));

%% Criterions

% Tangeante
Dt0phi=dt0phistreak(vect_t0,x0,t,vel,big_deltat0);
norm_Dt0phi = sqrt(sum(Dt0phi.^2,2));



% Tangeante correction
Dtiphi=dtiphistreak(vect_t0,x0,big_deltat0,t,vel);
norm_Dtiphi = sqrt(sum(Dtiphi.^2,2));

% % Time derivative of the tangeante norm
% Dt_norm_Dt0phi = ...
%     ( sqrt(sum( dt0phistreak(vect_t0(1:end-1),x0,t+deltat,vel,big_deltat0) .^2,2)) - ...
%     sqrt(sum( dt0phistreak(vect_t0(1:end-1),x0,t-deltat,vel,big_deltat0) .^2,2)) ) ...
%     / (2*deltat);
% Dt_norm_Dt0phi = Dt_norm_Dt0phi ./ norm_Dt0phi(1:end-1);
%
% % Time derivative of the tangeante correction norm
% Dt_norm_Dtiphi = ...
%     ( sqrt(sum( dtiphistreak(vect_t0(1:end-1),x0,big_deltat0,t+deltat,vel) .^2,2)) - ...
%     sqrt(sum( dtiphistreak(vect_t0(1:end-1),x0,big_deltat0,t-deltat,vel) .^2,2)) ) ...
%     / (2*deltat);
% Dt_norm_Dtiphi = Dt_norm_Dtiphi ./ norm_Dt0phi(1:end-1);
% Dt_norm_Dtiphi_Dtiphi = Dt_norm_Dtiphi ./ norm_Dtiphi(1:end-1);
%
% % Time derivative of the velocity norm (streamwise acceleration)
% % streakline
% streakline_mdt=streakline_2(x0,vect_t0(1:end-1),t-deltat,vel);
% streakline_pdt=streakline_2(x0,vect_t0(1:end-1),t+deltat,vel);
% % Norm of Eulerian velocity
% norm_v_mdt = sqrt(sum(vel(t-deltat,X,Y).^2,3));
% norm_v_pdt = sqrt(sum(vel(t+deltat,X,Y).^2,3));
% Dt_norm_V = ...
%     ( interpstreak_any(streakline_pdt,norm_v_pdt,X,Y) - ...
%     interpstreak_any(streakline_mdt,norm_v_mdt,X,Y) ) ...
%     / (2*deltat);
% % streakline_mdt=streakline(x0,t0,deltat0,t-deltat,vel);
% % streakline_pdt=streakline(x0,t0,deltat0,t+deltat,vel);
% % Dt_norm_V = ...
% %     ( sqrt(sum( interpstreak_any(streakline_pdt,v,X,Y,vect_t0) .^2,2)) - ...
% %       sqrt(sum( interpstreak_any(streakline_mdt,v,X,Y,vect_t0) .^2,2)) ) ...
% %       / (2*deltat);
% Dt_norm_V = Dt_norm_V ./ norm_Dt0phi(1:end-1);

% Normalized tangeante
n_Dt0phi = bsxfun(@times,1./norm_Dt0phi, Dt0phi);

% Lagrangian velocity
%norm_V = interpstreak_any(streakline1,norm_v,X,Y);
V = interpstreak_any(streakline1,v,X,Y);
norm_V = sqrt(sum(V.^2,2));

%% Remove points

iii = norm_Dt0phi < max(norm_Dt0phi(:))*2/3;

streakline1(iii,:)=[];
norm_Dt0phi(iii)=[];
norm_Dtiphi(iii)=[];
norm_V(iii)=[];

%% Guillermo idea

% % Meth 1 -> Equivalent to the tangant growth rate
% tic
% grad_tan_V = bsxfun(@times,1./(deltat0*norm_Dt0phi(1:end-1,1)), ...
%     (V(2:end,:)-V(1:end-1,:)) );
% strain_rate_on_streak = sum(n_Dt0phi(1:end-1,:) .* grad_tan_V,2);
% norm_grad_tan_V = sqrt(sum(grad_tan_V.^2,2));
% toc
% 
% % Meth 3
% tic
% grad_tan_V_3 = bsxfun(@times,1./(2*deltat0*norm_Dt0phi(2:end-1,1)), ...
%     (V(3:end,:)-V(1:end-2,:)) );
% strain_rate_on_streak_3 = sum(n_Dt0phi(2:end-1,:) .* grad_tan_V_3,2);
% norm_grad_tan_V_3 = sqrt(sum(grad_tan_V_3.^2,2));
% toc
% 
% % Meth2
% 
% tic
% nabla_v = gradient_mat_2(v,dX); % [ Mx My d n2 d]
% nabla_V = interpstreak_any(streakline1,nabla_v,X,Y); % [ N d n2 d]
% grad_tan_V_2 = squeeze(sum(bsxfun(@times,n_Dt0phi,nabla_V),2));
% strain_rate_on_streak_ref = sum(n_Dt0phi .* grad_tan_V_2,2);
% norm_grad_tan_V_2 = sum(grad_tan_V_2,2);
% toc

%% Pseudo Okubo-Weiss

% % Estimation of the vorticity
% n_Dt0phi_orto(:,1) = - n_Dt0phi(:,2);
% n_Dt0phi_orto(:,2) = n_Dt0phi(:,1);
% vort_estim = 2*sum(n_Dt0phi_orto(1:end-1,:) .* grad_tan_V,2);
% 
% % True vorticity
% omega=vort_mat(v,dX);%[Mx My d n2]
% omega = interpstreak_any(streakline1,omega,X,Y); % [ N d n2 d]
% 
% % Pseudo Okubo-Weiss
% pseudo_Q = ( strain_rate_on_streak.^2 - (2*vort_estim).^2 );


%%

ratio = 1;

% width=12;
% height=12;
% figure1=figure(1);
% set(figure1,'Units','inches', ...
%     'Position',[0 0 width height], ...
%     'PaperPositionMode','auto');
% 
% subplot(4,4,2)
% scatter(streakline1(:,1),streakline1(:,2),[],norm_Dt0phi,'filled')
% axis xy;axis equal;axis(ax);
% colorbar('location',loc_colorbar)
% cax=caxis;cax(2)=ratio*cax(2);
% caxis(cax);
% title('$ \| \partial_{t_0} \phi \| $',...
%     'FontUnits','points',...
%     'FontWeight','normal',...
%     'interpreter','latex',...
%     'FontSize',12,...
%     'FontName','Times');
% subplot(4,4,3)
% scatter(streakline1(:,1),streakline1(:,2),[],norm_Dtiphi,'filled')
% axis xy;axis equal;axis(ax);
% colorbar('location',loc_colorbar)
% caxis(cax);
% title('$\| \partial_{t_i} \phi \| $',...
%     'FontUnits','points',...
%     'FontWeight','normal',...
%     'interpreter','latex',...
%     'FontSize',12,...
%     'FontName','Times');
% subplot(4,4,4)
% scatter(streakline1(:,1),streakline1(:,2),[],norm_V,'filled')
% axis xy;axis equal;axis(ax);
% colorbar('location',loc_colorbar)
% caxis(cax);
% title('$\| V \| $',...
%     'FontUnits','points',...
%     'FontWeight','normal',...
%     'interpreter','latex',...
%     'FontSize',12,...
%     'FontName','Times');

%%

width=12;
height=12;
figure5=figure(5);
set(figure5,'Units','inches', ...
    'Position',[0 0 width height], ...
    'PaperPositionMode','auto');

% subplot(4,4,1)
% scatter(streakline1(:,1),streakline1(:,2),[],Qstreak,'filled')
% % plot(streakline1(:,1),streakline1(:,2),'o-')
% axis xy;axis equal;axis(ax);
% colorbar('location',loc_colorbar)
% title('$Q$',...
%     'FontUnits','points',...
%     'FontWeight','normal',...
%     'interpreter','latex',...
%     'FontSize',12,...
%     'FontName','Times');
% cax=caxis;
% subplot(4,4,2)
% scatter(streakline1(1:end-1,1),streakline1(1:end-1,2),[],pseudo_Q,'filled')
% % plot(streakline1(:,1),streakline1(:,2),'o-')
% axis xy;axis equal;axis(ax);
% caxis(cax);
% colorbar('location',loc_colorbar)
% title('pseudo-$Q$',...
%     'FontUnits','points',...
%     'FontWeight','normal',...
%     'interpreter','latex',...
%     'FontSize',12,...
%     'FontName','Times');
% 
% subplot(4,4,3)
% scatter(streakline1(1:end-1,1),streakline1(1:end-1,2),[],vort_estim,'filled')
% % plot(streakline1(:,1),streakline1(:,2),'o-')
% axis xy;axis equal;axis(ax);
% cax=caxis;
% colorbar('location',loc_colorbar)
% title('$\nabla^{\bot} \cdot v$ (estim)',...
%     'FontUnits','points',...
%     'FontWeight','normal',...
%     'interpreter','latex',...
%     'FontSize',12,...
%     'FontName','Times');
% subplot(4,4,4)
% scatter(streakline1(:,1),streakline1(:,2),[],omega,'filled')
% % plot(streakline1(:,1),streakline1(:,2),'o-')
% axis xy;axis equal;axis(ax);
% caxis(cax);
% colorbar('location',loc_colorbar)
% title('$\nabla^{\bot} \cdot v$ (ref)',...
%     'FontUnits','points',...
%     'FontWeight','normal',...
%     'interpreter','latex',...
%     'FontSize',12,...
%     'FontName','Times');



% subplot(4,4,5)
% scatter(streakline1(1:end-1,1),streakline1(1:end-1,2),[],...
%     strain_rate_on_streak,'filled')
% axis xy;axis equal;axis(ax);
% colorbar('location',loc_colorbar)
% title('$n^T S n $ (1st order)',...
%     'FontUnits','points',...
%     'FontWeight','normal',...
%     'interpreter','latex',...
%     'FontSize',12,...
%     'FontName','Times');
% cax=caxis;
% subplot(4,4,6)
% scatter(streakline1(2:end-1,1),streakline1(2:end-1,2),[],...
%     strain_rate_on_streak_3,'filled')
% axis xy;axis equal;axis(ax);
% colorbar('location',loc_colorbar)
% title('$n^T S n $ (2nd order)',...
%     'FontUnits','points',...
%     'FontWeight','normal',...
%     'interpreter','latex',...
%     'FontSize',12,...
%     'FontName','Times');
% caxis(cax);
% subplot(4,4,7)
% scatter(streakline1(:,1),streakline1(:,2),[],...
%     strain_rate_on_streak_ref,'filled')
% axis xy;axis equal;axis(ax);
% colorbar('location',loc_colorbar)
% caxis(cax);
% title('$n^T S n $ (ref)',...
%     'FontUnits','points',...
%     'FontWeight','normal',...
%     'interpreter','latex',...
%     'FontSize',12,...
%     'FontName','Times');


% subplot(4,4,9)
% scatter(streakline1(1:end-1,1),streakline1(1:end-1,2),[],...
%     Dt_norm_Dt0phi,'filled')
% axis xy;axis equal;axis(ax);
% colorbar('location',loc_colorbar)
% %cax=caxis;cax(2)=ratio*cax(2);
% cax=caxis;
% %caxis(cax);
% title('$\partial_t \|\partial_{t_0} \phi\| / \| \partial_{t_0} \phi \| $',...
%     'FontUnits','points',...
%     'FontWeight','normal',...
%     'interpreter','latex',...
%     'FontSize',12,...
%     'FontName','Times');
% subplot(4,4,10)
% scatter(streakline1(1:end-1,1),streakline1(1:end-1,2),[],...
%     Dt_norm_Dtiphi,'filled')
% axis xy;axis equal;axis(ax);
% colorbar('location',loc_colorbar)
% caxis(cax);
% title('$\partial_t \|\partial_{t_i} \phi \|/ \| \partial_{t_0} \phi \| $',...
%     'FontUnits','points',...
%     'FontWeight','normal',...
%     'interpreter','latex',...
%     'FontSize',12,...
%     'FontName','Times');
% subplot(4,4,11)
% scatter(streakline1(1:end-1,1),streakline1(1:end-1,2),[],...
%     Dt_norm_V,'filled')
% axis xy;axis equal;axis(ax);
% colorbar('location',loc_colorbar)
% caxis(cax);
% %title('Time derivative of the velocity norm')
% title('$\partial_t \|V \|/ \| \partial_{t_0} \phi \| $',...
%     'FontUnits','points',...
%     'FontWeight','normal',...
%     'interpreter','latex',...
%     'FontSize',12,...
%     'FontName','Times');
% subplot(4,4,12)
% scatter(streakline1(1:end-1,1),streakline1(1:end-1,2),[],...
%     Dt_norm_Dtiphi_Dtiphi,'filled')
% axis xy;axis equal;axis(ax);
% colorbar('location',loc_colorbar)
% caxis(cax);
% title('$\partial_t \|\partial_{t_i} \phi \|/ \| \partial_{t_i} \phi \| $',...
%     'FontUnits','points',...
%     'FontWeight','normal',...
%     'interpreter','latex',...
%     'FontSize',12,...
%     'FontName','Times');

subplot(3,1,1)
scatter(streakline1(:,1),streakline1(:,2),[],(norm_Dt0phi),'filled')
axis xy;axis equal;axis(ax);
colorbar('location',loc_colorbar)
cax=caxis;cax(2)=ratio*cax(2);
caxis(cax);
title('$ \| \partial_{t_0} \phi \| $',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'interpreter','latex',...
    'FontSize',12,...
    'FontName','Times');
subplot(3,1,2)
scatter(streakline1(:,1),streakline1(:,2),[],(norm_Dtiphi),'filled')
axis xy;axis equal;axis(ax);
colorbar('location',loc_colorbar)
caxis(cax);
title('$\| \partial_{t_i} \phi \| $',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'interpreter','latex',...
    'FontSize',12,...
    'FontName','Times');
subplot(3,1,3)
scatter(streakline1(:,1),streakline1(:,2),[],norm_V,'filled')
axis xy;axis equal;axis(ax);
colorbar('location',loc_colorbar)
caxis(cax);
title('$\| V \| $',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'interpreter','latex',...
    'FontSize',12,...
    'FontName','Times');

%%

% % figure(2);
% % plot(vect_t0,-V+Dtiphi,'r')
% % hold on;
% % plot(vect_t0,Dt0phi,'b')
% % hold off
% 
% % figure(3);
% % plot(vect_t0,grad_tan_V_2,'ro-')
% % hold on;
% % plot(vect_t0(2:end-1),grad_tan_V,'bo-')
% % hold off
% 
% % figure(3);
% % plot(vect_t0(1:end-1),Dt_norm_Dt0phi,'ro-')
% % hold on;
% % plot(vect_t0(1:end-1),strain_rate_on_streak,'bo-')
% % hold off
% 
% figure(3);
% plot(vect_t0(1:end-1),vort_estim,'ro-')
% hold on;
% plot(vect_t0,omega,'bo-')
% hold off
