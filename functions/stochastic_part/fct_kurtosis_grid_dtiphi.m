function [kurt,grid_streaks] = fct_kurtosis_grid_dtiphi...
    (grid,t0,deltat0,tau_adv,vel,grid_streaks,Dt0phi,ax)
%     (X,Y,t0,deltat0,tau_adv,vel)
% Compute the kurtosis of the streakline which is injected in (X,Y).
% The streakline is advected between t_0 and t0+tau_adv
% with a time step deltat0
%

%% Grid
x = grid.x_ref;
y = grid.y_ref;
[X,Y]=ndgrid(x,y);

%% Computation
t = t0 + tau_adv;

s=size(grid_streaks(:,:,:,1:end-1));
streakX=grid_streaks(:,:,1,1:end-1);
streakY=grid_streaks(:,:,2,1:end-1);
velocity_grid = vel(t,streakX(:),streakY(:));
velocity_grid = reshape(velocity_grid,[s(1:2) s(4) 2]);
velocity_grid = permute(velocity_grid,[1 2 4 3]);
Dtiphi = Dt0phi(:,:,:,1:end-1) + velocity_grid; clear Dt0phi

% clear grid_streaks;
Dtiphi(:,:,:,end+1) = 0;
norm_Dtiphi = sum(Dtiphi.^2,3); clear Dtiphi
kurt = squeeze(fct_kurtosis(permute(norm_Dtiphi,[4 1 2 3])));


% %% Computation
% t = t0 + tau_adv;
% grid_streaks=grid_of_streakline(X,Y,t0,deltat0,t,vel);
% Dt0phi = (grid_streaks(:,:,:,2:end)-grid_streaks(:,:,:,1:end-1))/deltat0;
% clear Dt0phi
% 
% s=size(grid_streaks(:,:,:,1:end-1));
% streakX=grid_streaks(:,:,1,1:end-1);
% streakY=grid_streaks(:,:,2,1:end-1);
% velocity_grid = vel(t,streakX(:),streakY(:));
% velocity_grid = reshape(velocity_grid,[s(1:2) s(4) 2]);
% velocity_grid = permute(velocity_grid,[1 2 4 3]);
% Dtiphi = Dt0phi + velocity_grid;
% 
% % clear grid_streaks;
% Dtiphi(:,:,:,end+1) = 0;
% norm_Dtiphi = sum(Dtiphi.^2,3); clear Dtiphi
% kurt = squeeze(fct_kurtosis(permute(norm_Dtiphi,[4 1 2 3])));

%% Pics on the sides
iii_remove = (X(:) == max(X(:)) | X(:) == min(X(:))) ...
    & (Y(:) == max(Y(:)) | Y(:) == min(Y(:)));
s=size(kurt);
kurt(iii_remove)=nan;
kurt=reshape(kurt,s);

%% Plots

%     ax=[0 2 0 1];

loc_colorbar = 'southoutside';
%colormap_ = 'default';
taille_police = 12;
width=6;
height=8;
% height=4;
% % width=12;
% % height=8;
X0=[0 0];
% X0=[8 0];
%X0=[5 20];
figure10=figure(10);
set(figure10,'Units','inches', ...
    'Position',[X0 width height], ...
    'PaperPositionMode','auto');
subplot(2,1,1)
imagesc(x,y,log(kurt'));
% imagesc(x,y,kurt'-3);
axis equal
axis xy
if nargin > 7
    axis(ax);
end
colorbar('location',loc_colorbar)
title('log(kurtosis$(\|\partial_{t_i} \phi\|)$)',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'interpreter','latex',...
    'FontSize',12,...
    'FontName','Times');
set(gca,...
    'Units','normalized',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',11,...
    'FontName','Times')
ylabel('x',...
    'FontUnits','points',...
    'interpreter','latex',...
    'FontSize',taille_police,...
    'FontName','Times')
xlabel('y',...
    'interpreter','latex',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',taille_police,...
    'FontName','Times')
%     axis(ax)
axis([x(1) x(end) y(1) y(end)])


subplot(2,1,2)
%imagesc(x,y,log(kurt'));
 imagesc(x,y,kurt'-3);
axis equal
axis xy
if nargin > 7
    axis(ax);
end
colorbar('location',loc_colorbar)
title('Excess kurtosis$(\|\partial_{t_i} \phi\|)$',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'interpreter','latex',...
    'FontSize',12,...
    'FontName','Times');
set(gca,...
    'Units','normalized',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',11,...
    'FontName','Times')
ylabel('x',...
    'FontUnits','points',...
    'interpreter','latex',...
    'FontSize',taille_police,...
    'FontName','Times')
xlabel('y',...
    'interpreter','latex',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',taille_police,...
    'FontName','Times')
%     axis(ax)
axis([x(1) x(end) y(1) y(end)])
