function  [X,Y] = phi_grid(X,Y,t0,tfinal,velocity,bool_plot,nsub,T_plot)
% Computation of the flow given an initial time and a grid of points(X,Y)
% t0 is the inital time, t the final time, velocity tje function which
% computes the fluid velocity at each time.
% bool_plot,nsub,T_plot are optional arguments for plots.
% If bool_plot = true, the functions plot the position of Lagrangian
% particles (spatially subsampled by a foactor nsub) and the vorticity 
% field every T_plot time steps
%

bool_mixing_diagno = true;

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
if bool_plot
    dX = [X(2,1)-X(1,1) Y(1,2)-Y(1,1)];
    X0=X;Y0=Y;
    ratio_increase = 1/30;
%     X_back = X; Y_back = Y;
    
    grid.dX=dX;
    grid.MX=size(X);
    grid.x_ref = X(:,1)';
    grid.y_ref = Y(1,:);
    
    folder_simu = [pwd '/images/' func2str(velocity) '/' ...
        'phi_grid_' ...
        '/Mx_' num2str(size(X0,1)) '_My_' num2str(size(X0,2)) '/'];
    mkdir(folder_simu);
    fprintf('Forward advection \n');
end

%% Choice of time step (CFL)
dt = fct_time_step(velocity);

%% Time of advection
N_t = ceil((tfinal-t0)/dt);

%% Time loop
for t=1:N_t
    % Time
    time = t0 + (t-1)*dt;
    
    % Forward advection (to get T(x_0,t) on a regular grid)
    [X,Y] = RK4_advection_lagrangienne(time, dt, X,Y, velocity);
    
%     if bool_plot
%         % Time
%         time_back = tfinal - (time -t0);
%         
%         % Forward advection (to get T(x_0,t) on a regular grid)
%         [X_back,Y_back] = RK4_advection_lagrangienne(time_back, - dt, ...
%             X_back,Y_back, velocity);
%     end
    
%     % Check if there are Nan values
%     if any(isnan(X(:)))
%         warning(['There are NaNs ' ...
%             '(possibly because some advected points left the domain).'...
%             'This behavior is expected for the wake flow.']);
%         %         error(['There are NaNs ' ...
%         %             '(possibly because some advected points left the domain)']);
%     end
    
    %% Plots
    if bool_plot && (mod(t,T_plot)==1)
        time
        t_string = num2str(time);
        t_string(t_string =='.')  = '_';
        
        Xsub = X(1:nsub:end,1:nsub:end);
        Ysub = Y(1:nsub:end,1:nsub:end);
        v = velocity(time,X0,Y0);
        vort = vort_mat(v,dX);
        figure(3)
        subplot(2,1,1)
        plot(Xsub(:),Ysub(:),'.');axis equal; axis xy; 
        if t==1
            axref=axis;
            switch func2str(velocity)
                case {'fct_wake_mega','fct_wake_megaRAM'}
                    axref(1)=0;                    
            end
            axref(1) = axref(1) - ratio_increase*(axref(2)-axref(1));
            axref(2) = axref(2) + ratio_increase*(axref(2)-axref(1));
            axref(3) = axref(3) - ratio_increase*(axref(4)-axref(3));
            axref(4) = axref(4) + ratio_increase*(axref(4)-axref(3));
        end
        axis(axref);
        title('Lagrangian particles');
        subplot(2,1,2)
        imagesc(X0(:,1),Y0(1,:),vort');
        axis xy;axis equal;colorbar;
        axis(axref);
        title('Vorticity')
        if t==1
            caxref=caxis;
            caxref=caxref/2;
            caxis(caxref);
        else
            caxis(caxref);
        end
        drawnow;pause(0.5);
        eval( ['print -depsc ' folder_simu ...
            'Lag_pcl_and_vort_' ...
            t_string '.eps']);
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


        %% Stretching and Mixing dignostic
        if bool_mixing_diagno
            % Gradient of the inverse (i.e. backward) flow
            nabla_phi = fct_nabla_phi(grid,X,Y);
            %nabla_phi = fct_nabla_phi(grid,X_back,Y_back);
            
            % Mezic criterion
            fct_mezic5(grid,nabla_phi,time - t0,axref);
            
            % FTLE
            fct_FTLE(grid,nabla_phi,time - t0,axref);
            
            % Okubo Weiss:
            % Estimation of the mean slpha^2 from the Okubo-Weiss assumptions
            % + plots in space
            fct_okubo_weiss(grid,v,axref);
            axis(axref);
        end
        
        %%
        drawnow;
    end
end
