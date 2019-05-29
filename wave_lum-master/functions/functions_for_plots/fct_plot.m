function [spectrum,name_plot,int_epsilon] = fct_plot...
    (model,fft_b_adv_part,day,ray_traject,ray_ampli)
% This function creates some plot online and save it
%

% Reduce buoyancy amplitude to obtain more realistic amplitudes
%(equivalent to reducing B0 and N by a factor 10)
model.odg_b = model.odg_b / 10;
fft_b_adv_part = fft_b_adv_part /10;

% idx_ = find(day=='_');
% day_title = eval(day(1:(idx_-1))) + eval(day((idx_+1):end))/100;
% day_title = num2str(day_title);
day_title = day;
day_title(day=='_') = '.';

%% Get paramters

% Grid
x = model.grid.x*1e-3;
y = model.grid.y*1e-3;
My = model.grid.MX(2);

% Other parameters
taille_police = 12;
id_part=1;
type_data = model.type_data;
folder_simu = model.folder.folder_simu;
plot_moments = model.advection.plot_moments;
%plot_epsilon_k = model.advection.plot_epsilon_k;
map = model.folder.colormap;

%% One particle
X0=[0 0];
T_adv_part = real(ifft2( fft_b_adv_part(:,:,1,id_part) ));

% if ( (eval(day) == 0) && ...
%         strcmp(model.type_data,'Perturbed_vortices') )
%     width = 3.2;
%     height = 3.2;
%     figure1=figure(1);
%     set(figure1,'Units','inches', ...
%         'Position',[X0(1) X0(2) width height], ...
%         'PaperPositionMode','auto');
%     contourf(x,y,T_adv_part');
%     x= model.grid.dX(1)*(0:model.grid.MX(1)-1)*1e-3;
%     y= model.grid.dX(2)*(0:model.grid.MX(2)-1)*1e-3;
%     Lx = model.grid.dX(1) * model.grid.MX(1);
%     sigma= 2 * Lx/15*1e-3;
%     center1x=x(1/4*model.grid.MX(1)+1);
%     center1y=y(1/4*model.grid.MX(2)+1);
%     nsig=40;
%     dist = 1.5;
%     rate = 0.3;
%     sigma = sigma/nsig;
%     center1x= 2e4*1e-3;
%     center1y=y(1/4*model.grid.MX(2)+1);
%     coord1=[center1x center1y];
%     size_square = 10e4*1e-3;
%     redline1 = [ [coord1(1)-size_square/2 coord1(2)-size_square/2] ; ...
%         [coord1(1)+size_square/2 coord1(2)-size_square/2] ; ...
%         [coord1(1)+size_square/2 coord1(2)+size_square/2] ; ...
%         [coord1(1)-size_square/2 coord1(2)+size_square/2] ; ...
%         [coord1(1)-size_square/2 coord1(2)-size_square/2] ];
%     hold on;
%     if strcmp(model.type_data,'Perturbed_vortices')
%         plot(redline1(:,1),redline1(:,2),'r','LineWidth',3);
%     elseif strcmp(model.type_data,'spot6')
%         plot(redline1(:,1),redline1(:,2),'b','LineWidth',3);
%     else
%         error('wrong type of data?');
%     end
%     center1x= x(1/2*model.grid.MX(1)+1) - 2e4*1e-3;
%     center1y=y(3/4*model.grid.MX(2)+1);
%     coord1=[center1x center1y];
%     redline1 = [ [coord1(1)-size_square/2 coord1(2)-size_square/2] ; ...
%         [coord1(1)+size_square/2 coord1(2)-size_square/2] ; ...
%         [coord1(1)+size_square/2 coord1(2)+size_square/2] ; ...
%         [coord1(1)-size_square/2 coord1(2)+size_square/2] ; ...
%         [coord1(1)-size_square/2 coord1(2)-size_square/2] ];
%     plot(redline1(:,1),redline1(:,2),'b','LineWidth',3);
%     hold off;
%
% else
width = 3.3;
height = 3.2;
figure1=figure(1);
close(figure1)
figure1=figure(1);
set(figure1,'Units','inches', ...
    'Position',[X0(1) 1 1.8*width 1.8*height], ...
    'PaperPositionMode','auto');
ax1 = axes;
[X,Y]=meshgrid(x,y);
s=surf(ax1,X,Y,T_adv_part);
s.EdgeColor = 'none';
view(2)
% end
axis(ax1,'xy')
axis(ax1 ,[x(1) x(end) y(1) y(end)])
caxis(ax1 ,[-1 1]*2*model.odg_b);
% caxis([-1 1]*1e-3);
set(ax1 ,...
    'Units','normalized',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',taille_police,...
    'FontName','Times')
ylabel(ax1 ,'y(km)',...
    'FontUnits','points',...
    'interpreter','latex',...
    'FontSize',taille_police,...
    'FontName','Times')
xlabel(ax1 ,'x(km)',...
    'interpreter','latex',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',taille_police,...
    'FontName','Times')
% title({'One realization', ...
%     ['\hspace{0.5cm} $t=' num2str(day) '$ day ']},...
title(ax1 ,['\hspace{0.5cm} $t=' num2str(day_title) '$ day '],...
    'FontUnits','points',...
    'FontWeight','normal',...
    'interpreter','latex',...
    'FontSize',12,...
    'FontName','Times')

% axis xy; axis equal
% % colormap(map)
% % colorbar


%% Lagrangian pcl
if model.advection.plot_charactericts
    %     ax1 = axes;
    
    %     hold on;
    x_ray_traject=ray_traject(:,1,:,:)*1e-3;
    y_ray_traject=ray_traject(:,2,:,:)*1e-3;
    % %     figure;
    %     plot(x_ray_traject(:),y_ray_traject(:),'.k', 'MarkerSize',20);
    %     hold on;
    ax2 = axes;
    plot(x_ray_traject(:),y_ray_traject(:),'.k', 'MarkerSize',20);
    hold on;
    plot(ax2,x_ray_traject(:),y_ray_traject(:),'.')
%     scatter(ax2,x_ray_traject(:),y_ray_traject(:),10,...
%         ones(size(y_ray_traject(:))),'filled')
%     scatter(ax2,x_ray_traject(:),y_ray_traject(:),10,ray_ampli(:),'filled')
    
    caxis(ax2 ,[0 4]);
    
    axis(ax2,'xy')
    axis(ax2,[x(1) x(end) y(1) y(end)])
    set(ax2 ,...
        'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontSize',taille_police,...
        'FontName','Times')
    ylabel(ax2 ,'y(km)',...
        'FontUnits','points',...
        'interpreter','latex',...
        'FontSize',taille_police,...
        'FontName','Times')
    xlabel(ax2 ,'x(km)',...
        'interpreter','latex',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontSize',taille_police,...
        'FontName','Times')
%     % title({'One realization', ...
%     %     ['\hspace{0.5cm} $t=' num2str(day) '$ day ']},...
%     title(ax2 ,['\hspace{0.5cm} $t=' num2str(day_title) '$ day '],...
%         'FontUnits','points',...
%         'FontWeight','normal',...
%         'interpreter','latex',...
%         'FontSize',12,...
%         'FontName','Times')
    
    %     axis xy; axis equal
    
    %%Link them together
    linkaxes([ax1,ax2])
    
    %%Give each one its own colormap
    colormap(ax1,map)
    colormap(ax2,'default')
    
    %%Hide the top axes
    ax2.Visible = 'off';
    ax2.XTick = [];
    ax2.YTick = [];
    %%Then add colorbars and get everything lined up
    %     axis xy;
    
    %
    %     plot(x_ray_traject(:),y_ray_traject(:),'.g', 'MarkerSize',1);
    % %     plot(x_ray_traject(:),y_ray_traject(:),'.k', 'MarkerSize',15);
    % %     plot(x_ray_traject(:),y_ray_traject(:),'.g');
    %     hold off;
    
    %%
    
    
    % %     axis xy;
    %     axis([ax1 ax2],'xy')
    %     axis([ax1 ax2],[x(1) x(end) y(1) y(end)])
    
    %     cb1 = colorbar(ax1);
    
    %     % caxis([-1 1]*1e-3);
    %     set(gca,...
    %         'Units','normalized',...
    %         'FontUnits','points',...
    %         'FontWeight','normal',...
    %         'FontSize',taille_police,...
    %         'FontName','Times')
    %     ylabel('y(km)',...
    %         'FontUnits','points',...
    %         'interpreter','latex',...
    %         'FontSize',taille_police,...
    %         'FontName','Times')
    %     xlabel('x(km)',...
    %         'interpreter','latex',...
    %         'FontUnits','points',...
    %         'FontWeight','normal',...
    %         'FontSize',taille_police,...
    %         'FontName','Times')
    %     % title({'One realization', ...
    %     %     ['\hspace{0.5cm} $t=' num2str(day) '$ day ']},...
    %     title(['\hspace{0.5cm} $t=' num2str(day_title) '$ day '],...
    %         'FontUnits','points',...
    %         'FontWeight','normal',...
    %         'interpreter','latex',...
    %         'FontSize',12,...
    %         'FontName','Times')
    % %     axis xy; axis equal
    %
    cb1 = colorbar(ax1,'Position',[.05 .11 .0675 .815]);
%     cb2 = colorbar(ax2,'Position',[.88 .11 .0675 .815]);
    
%     colorbar(ax1,'westoutside')
%     colorbar(ax2,'southoutside')
    
end

drawnow
eval( ['print -dpng ' folder_simu '/one_realization/'...
    num2str(day) '.png']);
% eval( ['print -depsc ' folder_simu '/one_realization/'...
%     num2str(day) '.eps']);

%% Spectrum
X0=[5 1];
% X0=[3.3 1];
close(figure(4))
figure4=figure(4);

widthtemp = 12;
heighttemp = 6;
set(figure4,'Units','inches', ...
    'Position',[X0(1) X0(2) widthtemp heighttemp], ...
    'PaperPositionMode','auto');
[spectrum,name_plot,int_epsilon] = fct_spectrum( model,fft_b_adv_part(:,:,:,id_part),'b');
set(gca,'XGrid','on','XTickMode','manual');
width = 4;
height = 3;
set(figure4,'Units','inches', ...
    'Position',[X0(1) X0(2) width height], ...
    'PaperPositionMode','auto');
set(gca,'YGrid','on')

set(gca,...
    'Units','normalized',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',taille_police,...
    'FontName','Times')
% ylabel('$|\hat{b}(\kappa)|^2$',...
ylabel('$E_b(\kappa) \bigl ( m^{3}.s^{-4}.{rad}^{-1} \bigr )$',...
    'FontUnits','points',...
    'interpreter','latex',...
    'FontSize',taille_police,...
    'FontName','Times')
xlabel('$\kappa \bigl ( rad.m^{-1} \bigr )$',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',taille_police,...
    'interpreter','latex',...
    'FontName','Times')
% title({'Spectrum of' ...
%     '\hspace{0.5cm} one realization'},...
title('Spectrum',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'interpreter','latex',...
    'FontSize',12,...
    'FontName','Times')
drawnow
eval( ['print -depsc ' folder_simu '/Spectrum/' day '.eps']);


%% Dissipation by scale
% if plot_epsilon_k
%     X0=[3.3 1];
%     close(figure(4))
%     figure4=figure(4);
%
%     widthtemp = 12;
%     heighttemp = 6;
%     set(figure4,'Units','inches', ...
%         'Position',[X0(1) X0(2) widthtemp heighttemp], ...
%         'PaperPositionMode','auto');
%     int_epsilon = fct_epsilon_k( model,fft_b_adv_part(:,:,:,id_part),...
%         int_epsilon_dt_m_1,'b');
%     set(gca,'XGrid','on','XTickMode','manual');
%     width = 4;
%     height = 3;
%     set(figure4,'Units','inches', ...
%         'Position',[X0(1) X0(2) width height], ...
%         'PaperPositionMode','auto');
%     set(gca,'YGrid','on')
%
%     set(gca,...
%         'Units','normalized',...
%         'FontUnits','points',...
%         'FontWeight','normal',...
%         'FontSize',taille_police,...
%         'FontName','Times')
%     ylabel('$\epsilon(\kappa)$',...
%         'FontUnits','points',...
%         'interpreter','latex',...
%         'FontSize',taille_police,...
%         'FontName','Times')
%     xlabel('$\kappa \bigl ( rad.m^{-1} \bigr )$',...
%         'FontUnits','points',...
%         'FontWeight','normal',...
%         'FontSize',taille_police,...
%         'interpreter','latex',...
%         'FontName','Times')
%     title({'Dissipation by scale' ...
%         '\hspace{0.5cm} one realization'},...
%         'FontUnits','points',...
%         'FontWeight','normal',...
%         'interpreter','latex',...
%         'FontSize',12,...
%         'FontName','Times')
%     drawnow
%     eval( ['print -depsc ' folder_simu '/Epsilon_k/' day '.eps']);
% end

%% Moments
if plot_moments
    T_adv_part = real(ifft2(fft_b_adv_part));
    tol = 1e-2;
    if strcmp(type_data,'Spectrum')
        tol = 1e-1;
    end
    mean_T = mean(T_adv_part,4);
    std_T = std(T_adv_part,0,4);
    odg_b = sqrt( mean(mean_T(:).^2) );
    
    % First and second order moments
    X0=[0 4.2];
    width = 3.65;
    height = 3;
    close(figure(2))%
    figure2=figure(2);
    set(figure2,'Units','inches', ...
        'Position',[X0(1) X0(2) 2*width height], ...
        'PaperPositionMode','auto');
    subplot(1,2,1)
    %     subimage(x,y,mean_T');axis xy;
    imagesc(x,y,mean_T');axis xy;
    axis equal
    caxis([-1 1]*model.odg_b);
    if model.folder.colormap_freeze
        colormap(map);
        colorbar;
        cbfreeze;
    else
        colormap(map);
        ax1 = gca;
        colorbar('peer',ax1);
        colormap(ax1,map);
    end
    set(gca,...
        'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontSize',taille_police,...
        'FontName','Times')
    ylabel('y(km)',...
        'FontUnits','points',...
        'interpreter','latex',...
        'FontSize',taille_police,...
        'FontName','Times')
    xlabel('x(km)',...
        'interpreter','latex',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontSize',taille_police,...
        'FontName','Times')
    title('\hspace{0.5cm} Mean',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'interpreter','latex',...
        'FontSize',12,...
        'FontName','Times')
    
    subplot(1,2,2)
    %     subimage(x,y,std_T');
    imagesc(x,y,std_T');axis xy;
    axis equal
    %caxis([0 model.odg_b/(1e-3)*1.5e-4]);
    if strcmp(type_data,'Spectrum')
        caxis([0 model.odg_b/(1e-3)*1e-3]);
    end
    if model.folder.colormap_freeze
        colormap('default');
        colorbar;
    else
        ax2 = gca;
        colorbar('peer',ax2);
        colormap(ax2,'default');
    end
    set(gca,...
        'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontSize',taille_police,...
        'FontName','Times')
    ylabel('y(km)',...
        'FontUnits','points',...
        'interpreter','latex',...
        'FontSize',taille_police,...
        'FontName','Times')
    xlabel('x(km)',...
        'interpreter','latex',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontSize',taille_police,...
        'FontName','Times')
    title('\hspace{1cm} Standard deviation',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'interpreter','latex',...
        'FontSize',12,...
        'FontName','Times')
    drawnow;
    eval( ['print -depsc ' folder_simu '/1st_2nd_order_moments/' day '.eps']);
    
    % Third and fourth order moments
    X0=[0 8.2];
    close(figure(3));%
    figure3=figure(3);
    set(figure3,'Units','inches', ...
        'Position',[X0(1) X0(2) 2*width height], ...
        'PaperPositionMode','auto');
    % These moments are shown only where the variance is large enough
    jjj = ( std_T < tol *odg_b );
    s=size(std_T);
    std_T=std_T(:);
    % The standard deviation used to compute skewness and kurtosis is set
    % to infinity when this standard deviation is too low
    % As such, skewness and kurtosis are equal to zero when this standard
    % deviation is too low
    std_T(jjj(:))=inf;
    std_T=reshape(std_T,s);
    
    % Centering
    T_prime = bsxfun(@plus,  real(ifft2( fft_b_adv_part)) , - mean_T);
    % Normalization
    T_prime = bsxfun(@times, T_prime, 1./ std_T);
    m3 = T_prime.^3 ;
    m3 = mean(m3,4);
    m4 = T_prime.^4 ;
    m4 = mean(m4,4);
    
    sm4=size(m4);
    m4=m4(:);
    m4(m4<3)=3;
    m4=reshape(m4,sm4);
    
    subplot(1,2,1)
    %     subimage(x,y,m3');
    imagesc(x,y,m3');axis xy;
    axis equal
    caxis(2*[-1 1]);
    if strcmp(type_data,'Spectrum')
        caxis(0.5*[-1 1]);
    end
    if model.folder.colormap_freeze
        colormap(map);
        colorbar;
        cbfreeze;
    else
        colormap(map);
        ax1 = gca;
        colorbar('peer',ax1);
        colormap(ax1,map);
    end
    set(gca,...
        'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontSize',taille_police,...
        'FontName','Times')
    ylabel('y(km)',...
        'FontUnits','points',...
        'interpreter','latex',...
        'FontSize',taille_police,...
        'FontName','Times')
    xlabel('x(km)',...
        'interpreter','latex',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontSize',taille_police,...
        'FontName','Times')
    title('\hspace{0.5cm} Skewness',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'interpreter','latex',...
        'FontSize',12,...
        'FontName','Times')
    
    subplot(1,2,2)
    %     subimage(x,y,log(m4'-3));
    imagesc(x,y,log(m4'-3));axis xy;
    axis equal
    if model.folder.colormap_freeze
        colormap('default');
        colorbar;
    else
        ax2 = gca;
        colorbar('peer',ax2);
        colormap(ax2,'default');
    end
    set(gca,...
        'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontSize',taille_police,...
        'FontName','Times')
    ylabel('y(km)',...
        'FontUnits','points',...
        'interpreter','latex',...
        'FontSize',taille_police,...
        'FontName','Times')
    xlabel('x(km)',...
        'interpreter','latex',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontSize',taille_police,...
        'FontName','Times')
    title('\hspace{1.5cm} log(Kurtosis-3)',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'interpreter','latex',...
        'FontSize',12,...
        'FontName','Times')
    if ~ strcmp(type_data, 'Spectrum')
        caxis([-5 3]);
    end
    drawnow
    eval( ['print -depsc ' folder_simu '/3rd_4th_order_moments/' day '.eps']);
    
end