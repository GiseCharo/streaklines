function  [X,Y] = phi_grid(X,Y,t0,t,velocity,bool_plot,nsub,T_plot)
% Computation of the flow given an initial time and a grid of points(X,Y)
% t0 is the inital time, t the final time, velocity tje function which
% computes the fluid velocity at each time.
% bool_plot,nsub,T_plot are optional arguments for plots.
% If bool_plot = true, the functions plot the position of Lagrangian
% particles (spatially subsampled by a foactor nsub) and the vorticity 
% field every T_plot time steps
%

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
end

%% Choice of time step (CFL)
dt = fct_time_step(velocity);

%% Time of advection
N_t = ceil((t-t0)/dt);

%% Time loop
for t=1:N_t
    % Time
    time = t0 + (t-1)*dt;
    
    % Forward advection (to get T(x_0,t) on a regular grid)
    [X,Y] = RK4_advection_lagrangienne(time, dt, X,Y, velocity);
    
    % Check if there are Nan values
    if any(isnan(X(:)))
        warning(['There are NaNs ' ...
            '(possibly because some advected points left the domain).'...
            'This behavior is expected for the wake flow.']);
        %         error(['There are NaNs ' ...
        %             '(possibly because some advected points left the domain)']);
    end
    
    %% Plots
    if bool_plot && (mod(t,T_plot)==1)
        Xsub = X(1:nsub:end,1:nsub:end);
        Ysub = Y(1:nsub:end,1:nsub:end);
        v = velocity(time,X0,Y0);
        vort = vort_mat(v,dX);
        figure(3)
        subplot(2,1,1)
        plot(Xsub(:),Ysub(:),'.');axis equal; axis xy; 
        if t==1
            axref=axis;
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
        else
            caxis(caxref);
        end
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
        drawnow;
    end
end
