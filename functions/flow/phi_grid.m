function  [X,Y] = phi_grid(X,Y,t0,t,velocity,bool_plot,nsub)
% Computation if the flow given an initial time and a regular grid of
% points
%

if nargin < 6
    bool_plot = false;
end
if nargin < 7
    nsub = 1;
end

%% Grid
% MX = size(X);
dX = [X(2,1)-X(1,1) Y(1,2)-Y(1,1)];

%% Choice of time step (CFL)
v = velocity(t0,X,Y);
dX=permute(dX,[1 3 2]);
dt=sum(bsxfun(@times,abs(v),pi*1./dX),3);
dt=max(dt(:));
dt=1/dt/2;
warning('This definition of the time step is not safe');

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
    
    %% PLots
    if bool_plot
        % nsub=1;
        Xsub = X(1:nsub:end,1:nsub:end);
        Ysub = Y(1:nsub:end,1:nsub:end);
        plot(Xsub(:),Ysub(:),'.');axis equal; axis xy; drawnow;
    end
end
