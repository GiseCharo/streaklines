clear;close all
%%computes the streakline
%%x0 is the point of injection,t0 the first time of injection and t the last,
%%deltat0 is the time interval of injection,


x0(1,1)=1;x0(1,2)=0.5;
% t0=0;deltat0=0.01;t=20;
t0=0;deltat0=1e-1;t_final=40;
% vel=@couetteplanS;
% vel=@DGyreS;
vel=@DGyreNS;
% vel=@couetteplanNS;
N= 50;
% x=linspace(0.75,1.25,2*N);
% y=linspace(0.25,0.75,N);
x=linspace(0,2,2*N);
y=linspace(0,1,N);

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


% % t0=0;deltat0=0.01;t=20;
% t0=0;deltat0=1e-1;t_final=10;
% N= 50;
% x=linspace(0,20,2*N);
% y=linspace(-6,6,N);
% vel=@fct_wake_megaRAM;

% Grid
grid.x_ref = x;
grid.y_ref = y;
% [X,Y]=ndgrid(x,y);

folder_res = ['images/' func2str(vel) '/jump_by_hand_everywhere/'];
mkdir([ folder_res 'dtiphi/'])
mkdir([ folder_res 'dt0phi/'])

for t=1:t_final
    t
    %% Computations
    %     grid_streaks=grid_of_streakline(X,Y,t0,deltat0,t,vel);
    %     Dt0phi = (grid_streaks(:,:,:,2:end)-grid_streaks(:,:,:,1:end-1))/deltat0;
    %     clear grid_streaks;
    %     Dt0phi(:,:,:,end+1) = - vel(t,X,Y);
    %     norm_Dt0phi = sum(Dt0phi.^2,3); clear Dt0phi
    %     kurt = squeeze(fct_kurtosis(permute(norm_Dt0phi,[4 1 2 3])));
    tic
    [~,grid_streaks,Dt0phi] = ...
        fct_jump_by_hand(grid,t0,deltat0,t-t0,vel);
    toc
    % kurt = fct_kurtosis_grid(X,Y,t0,deltat0,t-t0,vel);
    
    
    
    drawnow
    eval( ['print -depsc ' folder_res '/dt0phi/t=_' num2str(t) '.eps']);
    
    %         eval( ['print -depsc ' folder_simu ...
    %             'kurtosis_dt0phi_t0_' ...
    %             num2str(t0_local) '.eps']);
    
    % Kurtosis of the tangeante correction of the streakline
    tic
    fct_jump_by_hand_dtiphi...
        (grid,t0,deltat0,t-t0,vel,grid_streaks,Dt0phi);
    toc
    
    drawnow
    eval( ['print -depsc ' folder_res '/dtiphi/t=_' num2str(t) '.eps']);
    
    %         X = grid_streaks(:,:,1,1);
    %         Y = grid_streaks(:,:,2,1);
    
%     drawnow
%     %         eval( ['print -depsc ' folder_simu ...
%     %             'kurtosis_dtiphi_t0_' ...
%     %             num2str(t0_local) '.eps']);
    
    %% Plots
    
    % %     ax=[0 2 0 1];
    %
    %     loc_colorbar = 'southoutside';
    %     colormap_ = 'default';
    %     taille_police = 12;
    %     width=12;
    %     height=8;
    %     figure8=figure(8);
    %     set(figure8,'Units','inches', ...
    %         'Position',[5 20 width height], ...
    %         'PaperPositionMode','auto');
    %     imagesc(x,y,kurt'-3);
    %     axis equal
    %     axis xy
    %     colorbar('location',loc_colorbar)
    %     title('Excess kurtosis',...
    %         'FontUnits','points',...
    %         'FontWeight','normal',...
    %         'interpreter','latex',...
    %         'FontSize',12,...
    %         'FontName','Times');
    %     set(gca,...
    %         'Units','normalized',...
    %         'FontUnits','points',...
    %         'FontWeight','normal',...
    %         'FontSize',11,...
    %         'FontName','Times')
    %     ylabel('x',...
    %         'FontUnits','points',...
    %         'interpreter','latex',...
    %         'FontSize',taille_police,...
    %         'FontName','Times')
    %     xlabel('y',...
    %         'interpreter','latex',...
    %         'FontUnits','points',...
    %         'FontWeight','normal',...
    %         'FontSize',taille_police,...
    %         'FontName','Times')
    % %     axis(ax)
    %     axis([x(1) x(end) y(1) y(end)])
    
end