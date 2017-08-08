function [mean_nb_jump,grid_streaks,Dt0phi] = fct_jump_by_hand_dtiphi...
    (grid,t0,deltat0,tau_adv,vel,grid_streaks,Dt0phi,ax)
% (grid,t0,deltat0,tau_adv,vel,ax)
%     (X,Y,t0,deltat0,tau_adv,vel)
% Compute the kurtosis of the streakline which is injected in (X,Y).
% The streakline is advected between t_0 and t0+tau_adv
% with a time step deltat0
%


switch func2str(vel)
    case {'couetteplanS','couetteplanNS'}
        h=1;
        threshold = h/10;
        
    case {'DGyreS','DGyreNS'}
        L = 1;
        threshold = L/10;
        
    case {'fct_wake','fct_wake_megaRAM'}
        D = 1;
        threshold = D/10;
        
    otherwise
        error('Unknown function')
end
threshold = threshold/deltat0;

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
norm_Dtiphi = sqrt(sum(Dtiphi.^2,3)); clear Dtiphi

mean_nb_jump = mean(norm_Dtiphi > threshold,4);

%% Pics on the sides
iii_remove = (X(:) == max(X(:)) | X(:) == min(X(:))) ...
    & (Y(:) == max(Y(:)) | Y(:) == min(Y(:)));
s=size(mean_nb_jump);
mean_nb_jump(iii_remove)=nan;
mean_nb_jump=reshape(mean_nb_jump,s);

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
%X0=[5 20];
figure10=figure(10);
set(figure10,'Units','inches', ...
    'Position',[X0 width height], ...
    'PaperPositionMode','auto');
subplot(2,1,1)
imagesc(x,y,log(mean_nb_jump'));
% imagesc(x,y,kurt'-3);
axis equal
axis xy
if nargin > 7
    axis(ax);
end
colorbar('location',loc_colorbar)
title('log(Nb jump mean $(\|\partial_{t_i} \phi\|)$)',...
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
 imagesc(x,y,mean_nb_jump');
axis equal
axis xy
if nargin > 7
    axis(ax);
end
colorbar('location',loc_colorbar)
title('Nb jump mean $(\|\partial_{t_i} \phi\|)$',...
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
