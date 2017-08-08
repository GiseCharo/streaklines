function  mixing_criterions(X,Y,t0ini,tfinal,velocity,bool_plot,nsub,T_plot)
% Computation of the flow given an initial time and a grid of points(X,Y)
% t0ini is the inital time, t the final time, velocity tje function which
% computes the fluid velocity at each time.
% bool_plot,nsub,T_plot are optional arguments for plots.
% If bool_plot = true, the functions plot the position of Lagrangian
% particles (spatially subsampled by a foactor nsub) and the vorticity
% field every T_plot time steps
%

%% Parameters of the algortihm (to choose)
% delta t_0 / delta t
% ratio_dt0_dt = 1
ratio_dt0_dt = 10
%ratio_dt0_dt = 3

% tau / (pseudo-)period of the Eulerian velocity
ratio_tau_period = 4
%ratio_tau_period = 2
% ratio_tau_period = 0.5
%ratio_tau_period = 0.25
% ratio_tau_period = 1

% Which plot should appear
bool_back = false;
bool_plot_kurtosis = false;
bool_plot_OW = true;
bool_plot_Mezic = true;
bool_plot_vort = true;

bool_superimposed_pcl = true;
color_pcl = 'r';
MarkerSize = 20;

%% Plot parameters
if nargin < 6
    bool_plot = false;
end
if nargin < 7
    nsub = 1;
end
if nargin < 8
    T_plot = 10;
end
dX = [X(2,1)-X(1,1) Y(1,2)-Y(1,1)];
X0=X;Y0=Y;
ratio_increase = 1/30;

grid.dX=dX;
grid.MX=size(X);
grid.x_ref = X(:,1)';
grid.y_ref = Y(1,:);


%% Choice of time step (CFL)
dt = fct_time_step(velocity);

%% Choice of time step along the streaklines
% (no stability issues but a larger deltat0 implies a weaker convergence
% of the kurtosis)
deltat0 = ratio_dt0_dt * dt;

%% Times studied
N_t0 = ceil((tfinal-t0ini)/dt);

%% Time of advection
switch func2str(velocity)
    case {'fct_wake_megaRAM_2blocks','fct_wake_megaRAM','fct_wake',}
        tau_fixed = 5;
    case {'DGyreS','DGyreNS',}
        tau_fixed = 10;
    case {'couetteplanS','couetteplanNS'}
        h=1;
        nu=1;
        omega=2*nu*((2*pi/h)^2);
        tau_fixed = 2*pi/omega;
end
if bool_superimposed_pcl
    T_ini_pcl = tau_fixed;
    t0ini = t0ini + T_ini_pcl;
%     t0ini = t0ini + T_plot;
end

tau_fixed = ratio_tau_period * tau_fixed;
N_tau = ceil(tau_fixed/dt);


%%
X=X0;Y=Y0;
Xsub = X(1:nsub:end,1:nsub:end);
Ysub = Y(1:nsub:end,1:nsub:end);
Y0sub = Y0(1:nsub:end,1:nsub:end);
Xsub(abs(Y0sub)>0.5)=[];Ysub(abs(Y0sub)>0.5)=[];

v = velocity(t0ini,X0,Y0);
vort = vort_mat(v,dX);

figure(3)
subplot(2,1,1)
plot(Xsub(:),Ysub(:),'.');axis equal; axis xy;
% if t==1
axref=axis;
switch func2str(velocity)
    case {'fct_wake_mega','fct_wake_megaRAM'}
        axref(1)=0;
end
axref(1) = axref(1) - ratio_increase*(axref(2)-axref(1));
axref(2) = axref(2) + ratio_increase*(axref(2)-axref(1));
axref(3) = axref(3) - ratio_increase*(axref(4)-axref(3));
axref(4) = axref(4) + ratio_increase*(axref(4)-axref(3));
% end
axis(axref);
title('Lagrangian particles');
subplot(2,1,2)
imagesc(X0(:,1),Y0(1,:),vort');
axis xy;axis equal;colorbar;
axis(axref);
title('Vorticity')
% if t==1
caxref=caxis;
caxref = caxref/4;
caxis(caxref);
% else
%     caxis(caxref);
% end

folder_simu = [pwd '/images/' func2str(velocity) '/' ...
    'ratio_dt0_dt_' num2str(ratio_dt0_dt) ...
    'ratio_tau_period_' num2str(ratio_tau_period) ...
    '/Mx_' num2str(size(X0,1)) '_My_' num2str(size(X0,2)) '/'];
mkdir(folder_simu);


if bool_superimposed_pcl
    for t0_local_local=(t0ini-T_ini_pcl):dt:(t0ini-dt-T_plot)
        % Forward advection (to get T(x_0,t) on a regular grid)
        [Xsub,Ysub] = RK4_advection_lagrangienne(t0_local_local, ...
            dt, Xsub,Ysub, velocity);
    end
end

%% t_0 loop
for t0_idx=1:max([1 floor(T_plot/dt)]):N_t0
    t0_local = (t0_idx-1)*dt + t0ini;
    
    t0_local_string = num2str(t0_local);
    t0_local_string(t0_local_string =='.')  = '_';
    
    
    %% Advection along t_0
    if bool_superimposed_pcl
        for t0_local_local=(t0_local-T_plot):dt:(t0_local-dt)
            % Forward advection (to get T(x_0,t) on a regular grid)
            [Xsub,Ysub] = RK4_advection_lagrangienne(t0_local_local, ...
                dt, Xsub,Ysub, velocity);
        end
    end
    %%
    
    
    if ~bool_plot_kurtosis
        tfinal = t0_local + tau_fixed;
        
        % Initialisation
        X = X0; Y = Y0;
        if bool_back
            X_back = X0; Y_back = Y0;
        end
        
        %% tau loop
        for tau_local_idx=1:N_tau
            % Time
            tau_local = (tau_local_idx-1)*dt;
            time = t0_local + tau_local;
            
            % Forward advection (to get T(x_0,t) on a regular grid)
            [X,Y] = RK4_advection_lagrangienne(time, dt, X,Y, velocity);
            
            if bool_back
                % Time
                time_back = tfinal - tau_local;
                
                % Forward advection (to get T(x_0,t) on a regular grid)
                [X_back,Y_back] = RK4_advection_lagrangienne(time_back, - dt, ...
                    X_back,Y_back, velocity);
            end
            
        end
    end
    
    %% Kurtosis
    if bool_plot_kurtosis
        % Kurtosis of the tangeante of the streakline
        [~,grid_streaks,Dt0phi] = ...
            fct_kurtosis_grid(grid,t0_local,deltat0,tau_fixed,velocity,axref);
        %    [kurt,grid_streaks] = fct_kurtosis_grid(grid,X0,Y0,t0,deltat0,tau_adv,vel);
        
        X = grid_streaks(:,:,1,1);
        Y = grid_streaks(:,:,2,1);
        
        if bool_superimposed_pcl
            %             Xsub = X(1:nsub:end,1:nsub:end);
            %             Ysub = Y(1:nsub:end,1:nsub:end);
            %             Xsub(abs(Y0sub)>0.5)=[];Ysub(abs(Y0sub)>0.5)=[];
            subplot(2,1,1)
            hold on;plot(Xsub(:),Ysub(:),['.' color_pcl],'MarkerSize',MarkerSize);hold off;
            subplot(2,1,2)
            hold on;plot(Xsub(:),Ysub(:),['.' color_pcl],'MarkerSize',MarkerSize);hold off;
        end
        
        drawnow;pause(0.5);
        eval( ['print -depsc ' folder_simu ...
            'kurtosis_dt0phi_t0_' ...
            t0_local_string '.eps']);
        
        % Kurtosis of the tangeante correction of the streakline
        fct_kurtosis_grid_dtiphi...
            (grid,t0_local,deltat0,tau_fixed,velocity,grid_streaks,Dt0phi,axref);
        
        if bool_superimposed_pcl
            subplot(2,1,1)
            hold on;plot(Xsub(:),Ysub(:),['.' color_pcl],'MarkerSize',MarkerSize);hold off;
            subplot(2,1,2)
            hold on;plot(Xsub(:),Ysub(:),['.' color_pcl],'MarkerSize',MarkerSize);hold off;
        end
        
        drawnow;pause(0.5);
        eval( ['print -depsc ' folder_simu ...
            'kurtosis_dtiphi_t0_' ...
            t0_local_string '.eps']);
    end
    
    %% Lagrangian particles and vorticity
    %     if bool_plot && (mod(t,T_plot)==1)
    if bool_plot_vort
        %         Xsub = X(1:nsub:end,1:nsub:end);
        %         Ysub = Y(1:nsub:end,1:nsub:end);
        %         Xsub(abs(Y0sub)>0.5)=[];Ysub(abs(Y0sub)>0.5)=[];
        v = velocity(t0_local,X0,Y0);
        vort = vort_mat(v,dX);
        
        figure(3)
        subplot(2,1,1)
        plot(Xsub(:),Ysub(:),'.');axis equal; axis xy;
        %     if t==1
        %         axref=axis;
        %         switch func2str(velocity)
        %             case {'fct_wake_mega','fct_wake_megaRAM'}
        %                 axref(1)=0;
        %         end
        %         axref(1) = axref(1) - ratio_increase*(axref(2)-axref(1));
        %         axref(2) = axref(2) + ratio_increase*(axref(2)-axref(1));
        %         axref(3) = axref(3) - ratio_increase*(axref(4)-axref(3));
        %         axref(4) = axref(4) + ratio_increase*(axref(4)-axref(3));
        %     end
        axis(axref);
        title('Lagrangian particles');
        subplot(2,1,2)
        imagesc(X0(:,1),Y0(1,:),vort');
        axis xy;axis equal;colorbar;
        axis(axref);
        if bool_superimposed_pcl
            hold on;plot(Xsub(:),Ysub(:),['.' color_pcl],'MarkerSize',MarkerSize);hold off;
        end
        title('Vorticity')
        %     if t==1
        %         caxref=caxis;
        %     else
        caxis(caxref);
        %     end
        %         figure(4)
        %         subplot(2,1,1)
        %         imagesc(X0(:,1),Y0(1,:),v(:,:,1)');
        %         axis xy;axis equal;colorbar;
        %         title('Ux')
        %         axis(axref);
        %         subplot(2,1,2)
        %         imagesc(X0(:,1),Y0(1,:),v(:,:,1)');
        %         axis xy;axis equal;colorbar;
        %         title('Uy')
        %         axis(axref);
        
        drawnow;pause(0.5);
        eval( ['print -depsc ' folder_simu ...
            'vort_t0_' ...
            t0_local_string '.eps']);
    end
    
    
    %% Stretching and Mixing dignostic
    
    
    %     streakltemp = nan([size(X0) 2]);
    %     [X,Y]=phi_grid(X0,Y0,time0(i),t,velocity);
    %     streakltemp(:,:,1)=X;
    %     streakltemp(:,:,2)=Y;
    %     streakl(:,:,:,i)=streakltemp;
    
    % Gradient of the inverse (i.e. backward) flow
    nabla_phi = fct_nabla_phi(grid,X,Y);
    %nabla_phi = fct_nabla_phi(grid,X_back,Y_back);
    
    %% Mezic criterion
    if bool_plot_Mezic
        fct_mezic5(grid,nabla_phi,tau_fixed,axref);
        
        if bool_superimposed_pcl
            hold on;plot(Xsub(:),Ysub(:),['.' color_pcl],'MarkerSize',MarkerSize);hold off;
        end
        
        drawnow;pause(0.5);
        eval( ['print -depsc ' folder_simu ...
            'Mezic_t0_' ...
            t0_local_string '.eps']);
    end
    
    %% FTLE
    fct_FTLE(grid,nabla_phi,tau_fixed,axref);
    
    if bool_superimposed_pcl
        subplot(2,1,1)
        hold on;plot(Xsub(:),Ysub(:),['.' color_pcl],'MarkerSize',MarkerSize);hold off;
        subplot(2,1,2)
        hold on;plot(Xsub(:),Ysub(:),['.' color_pcl],'MarkerSize',MarkerSize);hold off;
    end
    
    drawnow;pause(0.5);
    eval( ['print -depsc ' folder_simu ...
        'FTLE_t0_' ...
        t0_local_string '.eps']);
    
    %% Okubo Weiss:
    if bool_plot_OW
        % Estimation of the mean slpha^2 from the Okubo-Weiss assumptions
        % + plots in space
        %v_back = velocity(time_back,X0,Y0);
        fct_okubo_weiss(grid,v,axref);
        axis(axref);
        
        %if bool_superimposed_pcl
        %hold on;plot(Xsub(:),Ysub(:),['.' color_pcl],'MarkerSize',MarkerSize);hold off;
        %end
        
        drawnow;pause(0.5);
        eval( ['print -depsc ' folder_simu ...
            'OW_t0_' ...
            t0_local_string '.eps']);
    end
    
    %%
    drawnow;
    %     end
end
