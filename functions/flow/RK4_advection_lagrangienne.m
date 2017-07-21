function [X,Y] = RK4_advection_lagrangienne(time, dt, X,Y, velocity)
% function [X,delta_X_per] = RK4_advection_lagrangienne(time, dt, X,Y, velocity)
% Inegration of the Lagrangian path X with a velcoity w
% using 4th order Runge-Kutta temporal scheme.
% The velocity w is defined by the function velocity
%

%% 4th order Runge-Kutta
k1 = velocity(time,X,Y);
k2 = velocity(time+dt/2,X+ k1(:,:,1)*dt/2,Y+ k1(:,:,2)*dt/2);
k3 = velocity(time+dt/2,X+ k2(:,:,1)*dt/2,Y+ k2(:,:,2)*dt/2);
k4 = velocity(time+dt,X+ k3(:,:,1)*dt,Y+ k3(:,:,2)*dt);

dX = (dt/3)*(k1/2 + k2 + k3 + k4/2);
X = X + dX(:,:,1);
Y = Y + dX(:,:,2);

% %% Deal with periodic boundaries conditions
% if model.advection.periodic_boundary_conditions
%     % Domain size
%     Lx=model.grid.dX(1)*model.grid.MX(1);
%     Ly=model.grid.dX(2)*model.grid.MX(2);
%     
%     % Save
%     X_true = X;
%     
%     % Keep X in the domain (in order to enable the evaluation of the
%     % velocity at this point in the next iteration)
%     X(:,1)= mod(X(:,1),Lx);
%     X(:,2)= mod(X(:,2),Ly);
%     
%     % Save the difference
%     % (in order to be able to compute the flow gradient)
%     delta_X_per = X_true - X;
%     clear X_true
% else
%     delta_X_per = 0;
% end
