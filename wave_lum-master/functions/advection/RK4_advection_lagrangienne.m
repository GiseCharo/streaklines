function [wave_group,d_wave_group,delta_wave_group_per] = ...
    RK4_advection_lagrangienne(advection,grid, ...
    wave_group, ratio_ampli_ini, w_Eul, grad_w_Eul, w_Lag, X0)
% function [X,dX,delta_X_per] = RK4_advection_lagrangienne(advection,grid, X, w, X0)
% Inegration of the Lagrangian path X with a velcoity w
% using 4th order Runge-Kutta temporal scheme.
% The velocity w is defined on the grid X0
%

%% Time step
dt=advection.dt_adv;

%% Replicate on borders to enable the interpolation on characteristics
if advection.periodic_boundary_conditions
    w_Eul = fct_mirror_on_border_period(w_Eul,advection.nbp);
    grad_w_Eul = fct_mirror_on_border_period(grad_w_Eul,advection.nbp);
else
    w_Eul = fct_mirror_on_border(w_Eul,advection.nbp);
    grad_w_Eul = fct_mirror_on_border(grad_w_Eul,advection.nbp);
end

%% Time integration
k1 = deriv_advection_lagrangienne(advection, wave_group, ...
    w_Eul, grad_w_Eul, w_Lag, X0);
if advection.sto
    % Euler–Maruyama
    d_wave_group = dt*k1;
else
    % 4th order Runge-Kutta
    k2 = deriv_advection_lagrangienne(advection, wave_group + k1*dt/2, ...
        w_Eul, grad_w_Eul, w_Lag, X0);
    k3 = deriv_advection_lagrangienne(advection, wave_group + k2*dt/2, ...
        w_Eul, grad_w_Eul, w_Lag, X0);
    k4 = deriv_advection_lagrangienne(advection, wave_group + k3*dt, ...
        w_Eul, grad_w_Eul, w_Lag, X0);
    
    d_wave_group = (dt/3)*(k1/2 + k2 + k3 + k4/2);
end
% warning('the wave vector time integration must be done inside here');
wave_group = wave_group + d_wave_group;

%% Update Amplitude
% wave_group(:,5,:,:) = ratio_ampli_ini ...
%    .* ( wave_group(:,3,:,:).^2 + wave_group(:,4,:,:).^2 ).^(1/4) ;

%% Deal with periodic boundaries conditions
if advection.periodic_boundary_conditions
    % Domain size
    Lx=grid.dX(1)*grid.MX(1);
    Ly=grid.dX(2)*grid.MX(2);
    
    % Save
    wave_group_true = wave_group(:,1:2,:,:);
    
    % Keep X in the domain (in order to enable the evaluation of the
    % velocity at this point in the next iteration)
    wave_group(:,1,:,:)= mod(wave_group(:,1,:,:),Lx);
    wave_group(:,2,:,:)= mod(wave_group(:,2,:,:),Ly);
    
    % Save the difference
    % (in order to be able to compute the flow gradient)
    delta_wave_group_per = wave_group_true - wave_group(:,1:2,:,:);
    clear X_true
else
    delta_wave_group_per = 0;
end

%% Sub function
    function dwave_group_=deriv_advection_lagrangienne...
            (advection_, wave_group_,w_Eul, grad_w_Eul, w_Lag, X0_)
        % Interpolate the velcoity on the Lagrangian particles
        %
        
        % Initialization
        siz = size(wave_group_);
        dwave_group_ = nan(siz);
        grad_w_Eul_on_wave = nan([siz(1) 2 2 advection_.N_ech]);
        
        if size(w_Eul,4)>1 % Random Eulerian velocity field
            if size(w_Eul,4) ~= advection_.N_ech
                error(['The number of velocity realization does not match' ...
                    ' the number of ray realizations']);
            end
            for k=1:advection_.N_ech
                % Add eulerian velocity interpolation
                for p=1:2
                    dwave_group_(:,p,:,k) = ...
                        interp2(X0_{1},X0_{2},w_Eul(:,:,p,k)',...
                        wave_group_(:,1,:,k),wave_group_(:,2,:,k));
                end
                % Eulerian velocity gradient interpolation
                for p=1:2
                    for q=1:2
                        grad_w_Eul_on_wave(:,p,q,k) = ...
                            interp2(X0_{1},X0_{2},grad_w_Eul(:,:,p,k,q)',...
                            wave_group_(:,1,:,k),wave_group_(:,2,:,k));
                    end
                end
            end
        else % Deterministic Eulerian velocity field
            for k=1:advection_.N_ech
                % Add eulerian velocity interpolation
                for p=1:2
                    dwave_group_(:,p,:,k) = ...
                        interp2(X0_{1},X0_{2},w_Eul(:,:,p)',...
                        wave_group_(:,1,:,k),wave_group_(:,2,:,k));
                end
                % Eulerian velocity gradient interpolation
                for p=1:2
                    for q=1:2
                        grad_w_Eul_on_wave(:,p,q,k) = ...
                            interp2(X0_{1},X0_{2},grad_w_Eul(:,:,p,:,q)',...
                            wave_group_(:,1,:,k),wave_group_(:,2,:,k));
                    end
                end
            end
        end
        % Add Lagrangian velocity interpolation
        dwave_group_(:,1:2,:,:) = dwave_group_(:,1:2,:,:) + w_Lag;
        
%         % Wave vector dynamics
%         dwave_group_(:,3:4,:,:) = ...
%             - grad_w_Eul_on_wave(:,:,1,:).* wave_group_(:,3,:,:) ...
%             - grad_w_Eul_on_wave(:,:,2,:).* wave_group_(:,4,:,:) ;
%         clear grad_w_Eul_on_wave
%         
%         % Wave amplitude dynamics
%         % ("slaved" update outside the function)
%         dwave_group_(:,5,:,:) = 0;
    end
end
