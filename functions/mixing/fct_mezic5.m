function fct_mezic5(grid,nabla_meso_v,time_t,ax)
% function fct_mezic5(grid,nabla_meso_v,time_t,dot_red)
% This function plots Mezic's mixing criterion as well as personnal
% carcterisation of evolution of gradient
%

loc_colorbar = 'southoutside';

%% From the gradient of the flow to the  gradient of mesochronic velocity

% Remove the identity matrix
nabla_meso_v(:,:,1,1) = nabla_meso_v(:,:,1,1)-1;
nabla_meso_v(:,:,2,2) = nabla_meso_v(:,:,2,2)-1;

% Divide by the time of advection
nabla_meso_v = 1/ time_t * nabla_meso_v;

% Trace of the mesochronic velocity gradient
tr_meso_v = nabla_meso_v(:,:,1,1) + nabla_meso_v(:,:,2,2);

% Determinant of the mesochronic velocity gradient
det_meso_v = - 1/ time_t * tr_meso_v;

%% Mixing criterion
criterion1 = det_meso_v < 0; % Mesohyperbolic without flipping
criterion2 = (det_meso_v >= 0) & (det_meso_v <= 4/time_t^2); % Mesoelliptic
criterion3 = det_meso_v > 4/time_t^2; % Mesohyperbolic with flipping
criterion = 1.8*criterion1 + 2.5 * criterion2 + 3 * criterion3;
% If the two types of mesohyperbolicity (criterion =
% 1.8 and criterion = 3) are next to each other, there is mixing.

%% Plot

% day = num2str(floor(t*model.advection.dt_adv/(3600*24)));

% Grid
x = grid.x_ref;
y = grid.y_ref;

width=4;

%X0=[8 3];
 X0=[4 3];
height=4;
figure9=figure(9);
set(figure9,'Units','inches', ...
    'Position',[X0 width height], ...
    'PaperPositionMode','auto');
imagesc(x,y,criterion')
axis xy
axis equal
title('Mezic states',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'interpreter','latex',...
    'FontSize',12,...
    'FontName','Times');
set(gca,...
    'Clim',[0 4], ...
    'Units','normalized',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',11,...
    'FontName','Times')
ylabel('y',...
    'FontUnits','points',...
    'interpreter','latex',...
    'FontSize',11,...
    'FontName','Times')
xlabel('x',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',11,...
    'FontName','Times')
axis([x(1) x(end) y(1) y(end)])
if nargin > 3
    axis(ax);
%     hold on;
%     plot(dot_red(1,:),dot_red(2,:),'.r');
end
colorbar('location',loc_colorbar)
%drawnow
% eval( ['print -depsc ' model.folder.folder_simu '/mezic_state/' ...
%     num2str(day) '.eps']);


