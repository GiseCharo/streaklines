function  [X_back,Y_back] = phi_grid_back(X_back,Y_back,t0,tfinal,velocity,bool_plot,nsub,T_plot)
% Computation of the backward flow given an final time and a grid of points(X_back,Y_back)
% t0 is the inital time, tfinal the final time, velocity the function which
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
    dX = [X_back(2,1)-X_back(1,1) Y_back(1,2)-Y_back(1,1)];
    X0=X_back;Y0=Y_back;
    ratio_increase = 1/30;
    
    grid.dX=dX;
    grid.MX=size(X_back);
    grid.x_ref = X_back(:,1)';
    grid.y_ref = Y_back(1,:);
    
    folder_simu = [pwd '/images/' func2str(velocity) '/' ...
        'phi_grid_' ...
        '/Mx_' num2str(size(X0,1)) '_My_' num2str(size(X0,2)) '/'];
    mkdir(folder_simu);
end

%% Choice of time step (CFL)
dt = fct_time_step(velocity);

%% Time of advection
N_t = ceil((tfinal-t0)/dt);
fprintf('Backward advection \n');

%% Time loop
for t=1:N_t
    % Time
    time = t0 + (t-1)*dt;
    
%     % Forward advection (to get T(x_0,t) on a regular grid)
%     [X,Y] = RK4_advection_lagrangienne(time, dt, X,Y, velocity);
%
%     if bool_plot
% Time
time_back = tfinal - (time -t0);

% Forward advection (to get T(x_0,t) on a regular grid)
[X_back,Y_back] = RK4_advection_lagrangienne(time_back, - dt, ...
    X_back,Y_back, velocity);
%     end
    
%     % Check if there are Nan values
%     if any(isnan(X_back(:)))
%         warning(['There are NaNs ' ...
%             '(possibly because some advected points left the domain).'...
%             'This behavior is expected for the wake flow.']);
%         %         error(['There are NaNs ' ...
%         %             '(possibly because some advected points left the domain)']);
%     end
    
    %% Plots
    if bool_plot && (mod(t,T_plot)==1)
        time_back
        t_string = num2str(time_back);
        t_string(t_string =='.')  = '_';
        
        Xsub = X_back(1:nsub:end,1:nsub:end);
        Ysub = Y_back(1:nsub:end,1:nsub:end);
        v = velocity(time_back,X0,Y0);
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
            % nabla_phi = fct_nabla_phi(grid,X,Y);
            nabla_phi = fct_nabla_phi(grid,X_back,Y_back);
            
            % Mezic criterion
            fct_mezic5(grid,nabla_phi,tfinal-time_back,axref);
            
            % FTLE
            fct_FTLE(grid,nabla_phi,tfinal-time_back,axref);
            
            % Okubo Weiss:
            % Estimation of the mean alpha^2 from the Okubo-Weiss assumptions
            % + plots in space
            v_back = velocity(time_back,X0,Y0);
            fct_okubo_weiss(grid,v_back,axref);
            axis(axref);
        end
        
        %%
        drawnow;
    end
end
