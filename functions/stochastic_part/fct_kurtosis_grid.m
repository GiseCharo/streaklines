function [kurt,grid_streaks,Dt0phi] = fct_kurtosis_grid...
    (grid,t0,deltat0,tau_adv,vel,ax)
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
grid_streaks=grid_of_streakline(X,Y,t0,deltat0,t,vel);
Dt0phi = (grid_streaks(:,:,:,2:end)-grid_streaks(:,:,:,1:end-1))/deltat0;
% clear grid_streaks;
Dt0phi(:,:,:,end+1) = - vel(t,X,Y);
norm_Dt0phi = sum(Dt0phi.^2,3); 
%clear Dt0phi
kurt = squeeze(fct_kurtosis(permute(norm_Dt0phi,[4 1 2 3])));

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
X0=[8 0];
%X0=[5 20];
figure8=figure(8);
set(figure8,'Units','inches', ...
    'Position',[X0 width height], ...
    'PaperPositionMode','auto');
subplot(2,1,1)
imagesc(x,y,log(kurt'));
% imagesc(x,y,kurt'-3);
axis equal
axis xy
if nargin > 5
    axis(ax);
end
colorbar('location',loc_colorbar)
title('log(kurtosis$(\|\partial_{t_0} \phi\|)$)',...
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
if nargin > 5
    axis(ax);
end
colorbar('location',loc_colorbar)
title('Excess kurtosis$(\|\partial_{t_0} \phi\|)$',...
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
