function [fft_b, model] = fct_fft_advection_sto_mat(model,  fft_b)
% Advection of buoyancy using SQG or SQG_MU model
%

tic
%% Folder to save plots and files
init_folder_to_save

%% Grid

% Spatial grid
model.grid.x = model.grid.dX(1)*(0:model.grid.MX(1)-1);
model.grid.y = model.grid.dX(2)*(0:model.grid.MX(2)-1);

% Grid in Fourier space
model = init_grid_k (model);

% Ensemble size
% N_ech=model.advection.N_ech;

%% Preprocessing of the Lagrangian particles

% For now, we do not consider masks since we assume periodic boudary
% conditions

% Reference grid
X0{1} = model.grid.dX(1)* ...
    ((-model.wave.advection.nbp):...
      (model.grid.MX(1)+model.wave.advection.nbp-1));
X0{2} = model.grid.dX(2)* ...
    ((-model.wave.advection.nbp):...
      (model.grid.MX(2)+model.wave.advection.nbp-1));
% X0{1} = repmat( model.grid.dX(1)* ((-nbp):(MX(1)+nbp-1)) ...
%     ,[1 1 1 model.advection.N_ech]);
% X0{2} = repmat( model.grid.dX(2)* ((-nbp):(MX(2)+nbp-1)) ...
%     ,[1 1 1 model.advection.N_ech]);

wave_group_ini = init_Lagrangian_pcl(model);
wave_group = wave_group_ini;
wave_group_forward = wave_group_ini;
% ratio_ampli_ini = wave_group_ini(:,5,:,:) ...
%    .* ( wave_group_ini(:,3,:,:).^2 + wave_group_ini(:,4,:,:).^2 ).^(-1/4) ;

if model.advection.plot_charactericts
    ray_wave_group = wave_group;
end

%% Initialisation of the spatial fields

% Create several identical realizations of the intial buoyancy
fft_b = repmat(fft_b(:,:,1),[1 1 1 model.advection.N_ech]);

% Remove aliasing
fft_b(model.grid.k.ZM(1),:,:,:)=0;
fft_b(:,model.grid.k.ZM(2),:,:)=0;

% Initial large-scale velocity
fft_w = SQG_large_UQ(model, fft_b);
w=real(ifft2(fft_w));

% figure;imagesc(real(ifft2(fft_b))');axis xy;axis equal;colorbar;
% figure;imagesc(w(:,:,1)');axis xy;axis equal;colorbar;

% % Create several identical realizations of the intial buoyancy
% fft_b = repmat(fft_b(:,:,1),[1 1 1 model.advection.N_ech]);

%% Choice of the variance tensor a
model = fct_init_sigma(model,fft_w);

%% Forcing
model = fct_set_forcing(model);

%% Hyperviscosity
model = fct_set_HV(model,w);

%% Choice of time step : CFL
model.advection.dt_adv = fct_CFL(model.grid,model.advection,model.sigma,w);
% model.advection.dt_adv = fct_CFL(model,w);

%% Statistics of the time-uncorrelated velocity
[model,sigma_adv,coef_modulation,sigma_s,sigma_n_s] = ...
    fct_sigma(model,fft_w);

%% Loop on time

% Used for the first plot
tt_last = -inf;

% Number of time step
% N_t = ceil(model.advection.advection_duration/model.advection.dt_adv);

first_wave_plot = true;
fft_wave_noise = nan;

%% Print informations

if model.plots
    % Printing some information
    fprintf(['The initial condition is ' model.type_data ' \n'])
    fprintf(['1/k_c is equal to ' num2str(1/model.sigma.k_c) ' m \n'])
    fprintf(['Time step : ' num2str(model.advection.dt_adv) ' seconds \n']);
    fprintf(['Time of advection : ' num2str(...
        model.advection.advection_duration/3600/24) ' days \n']);
    fprintf(['Ensemble size : ' num2str(model.sigma.N_ech) ' realizations \n']);
    fprintf(['Resolution : ' num2str(model.grid.MX(1)) ' x ' ...
        num2str(model.grid.MX(2)) ' \n']);
    if model.advection.HV.bool
        str_subgridtensor = 'Hyper-viscosity';
    elseif model.advection.Lap_visco.bool
        str_subgridtensor = 'Laplacian diffusion';
    else
        str_subgridtensor = 'None';
    end
    if ( model.advection.HV.bool | model.advection.Lap_visco.bool) & ...
            model.advection.Smag.bool
        str_subgridtensor = ['Heterogneous ' str_subgridtensor];
    end
    fprintf(['Deterministic subgrid tensor : ' str_subgridtensor ' \n']);
    fprintf(['Model type : ' add_subgrid_deter ' \n']);
    if model.sigma.sto | model.advection.Smag.bool
        fprintf(['Details : ' subgrid_details ' \n']);
    end
    if model.sigma.sto
        fprintf(['type spectrum sigma :' model.sigma.type_spectrum ' \n']);
    end
end

if model.sigma.sto & ...
        strcmp(model.sigma.type_spectrum,'SelfSim_from_LS') & ...
        model.sigma.time_smooth.bool
    sigma_s = 0;
    % sigma_s = sigma;
end

%% Use a saved files of a former simulation ?
if model.advection.use_save
    warning(['The run begin from an older file instead of from the' ...
        'initial condition']);
    day = num2str(model.advection.day_save);
    name_file = [model.folder.folder_simu '/files/' day '.mat']; clear day
    clear fft_b
    model_ref = model;
    
    load(name_file)
    %      model.advection.advection_duration =  ...
    %          model_ref.advection.advection_duration;
    %      model = catstruct(model_ref,model);
    model = model_ref;
    
    if (exist('fft_b','var')==1)
    elseif (exist('fft_buoy_part','var')==1)
        fft_b = fft_buoy_part; clear fft_buoy_part
    elseif (exist('fft_T_adv_part','var')==1)
        fft_b = fft_T_adv_part; clear fft_T_adv_part
    elseif (exist('fft_T_adv_part','var')==1)
        fft_b = fft_T_adv_part; clear fft_T_adv_part
    elseif (exist('fft_tracer_part','var')==1)
        fft_b = fft_tracer_part; clear fft_tracer_part
    elseif (exist('fft_buoy_part_ref','var')==1)
        fft_b = fft_buoy_part_ref; clear fft_buoy_part_ref
    else
        error('Cannot find buoyancy field')
    end
    if model.sigma.sto
%         if size(fft_b,4) < model.sigma.N_ech
        if size(fft_b,4) < model.advection.N_ech
            if size(fft_b,4) == 1
                fft_b = repmat ( fft_b, [ 1 1 1 model.advection.N_ech ]);
                clear w fft_w
            else
                error('The number of realisation of the saved file is too low');
            end
        end
        if size(fft_b,4) > model.advection.N_ech
            warning(['The number of realisation of the saved file is too high.' ...
                ' Some realisations are hence removed.']);
            fft_b(:,:,:,(model.advection.N_ech+1):end) = [];
        end
        if ( strcmp(model.sigma.type_spectrum,'EOF') || ...
                strcmp(model.sigma.type_spectrum,'Euler_EOF') ) ...
                && ( model.sigma.nb_EOF > model.advection.N_ech )
            warning(['The number of EOF is larger than the ensemble size.' ...
                ' Some EOFs are hence removed.']);
            model.sigma.nb_EOF = model.advection.N_ech;
            sigma(:,:,:,:,(model.advection.N_ech+1):end) = [];
        end
    end
    
    % Version of matlab
    vers = version;
    year = str2double(vers(end-5:end-2));
    subvers = vers(end-1);
    model.folder.colormap_freeze = ...
        (  year < 2014 || ( year == 2014 && strcmp(subvers,'a') ) );
    
    if ~isfield(model,'plots')
        model.plots = true;
    end
    
    %     t_ini=t+1;
    %     time = time
    
    if ~(exist('time','var')==1)
        time =t*model.advection.dt_adv;
    end
    
    day_last_plot = floor(time/24/3600);
else
    %     t_ini=1;
    
    
    if model.advection.cov_and_abs_diff
        cov_w = nan(1,N_t);
        day_ref_cov = 100;
        %day_ref_cov = 1;
        %warning(['This reference time for the covariance should ' ...
        %    'correspond to the stationnary regime']);
        % w_save = w;
        load( [model.folder.folder_simu '/files/' num2str(day_ref_cov)...
            '.mat'],'w','t');
        w_ref_cov = w(:)'; clear w;
        t_ref_cov = t; clear t;
        % w = w_save; clear w_save;
    end
    
    time = 0;
    w_fv = w;
    
    day_last_plot = -inf;
end

% if model.sigma.sto
%     for sampl=1:N_ech
%         model_sampl(sampl)=model;
%     end
% end

%%

bool_save = true;

while time < model.advection.advection_duration
    %% Time-correlated velocity
    fft_w = SQG_large_UQ(model, fft_b);
    w = real(ifft2(fft_w));
    
    %% Time-uncorrelated velocity 
    if  model.advection.sto
        [model,sigma_adv,coef_modulation,sigma_s,sigma_n_s] = ...
            fct_sigma(model,fft_w);
    end
%     sigma_dBt_on_sq_dt = ...
%         fct_sample_sigma_dBt_on_sq_dt(model,sigma_adv,coef_modulation);
    
    %% Specify determinstic heterogeneous subgrid model
    if model.advection.Smag.bool
        if model.advection.Lap_visco.bool
            % Heterogeneous dissipation coefficient
            model.advection.coef_diff = fct_coef_diff(model,fft_b);
            % Coefficient coef_Smag to target a specific diffusive scale
            model.advection.coef_diff = ...
                model.advection.Smag.coef_Smag * ...
                model.advection.coef_diff ;
        elseif model.advection.HV.bool
            % Heterogeneous HV coefficient
            model.advection.coef_diff = ...
                fct_coef_HV(model,fft_b);
        end
        % Maximum dissipation coefficient
        model.advection.HV.maxVal = max(model.advection.coef_diff(:));
    end
    
    %% Dynamical time step
    model.advection.dt_adv = fct_CFL(model.grid,model.advection,...
        model.sigma,w);
    
    %% Ray integration Forward
    
    % Definition of currents for the wave dynamics
    fft_for_wave = fft_w;
    if model.wave.advection.smooth_LS_current
        % Filtering
        fft_for_wave((model.wave.grid.MX(1)/2+1): ...
            (end-1-model.wave.grid.MX(1)/2+1),:,:,:)=0;
        fft_for_wave(:,(model.wave.grid.MX(2)/2+1): ...
            (end-1-model.wave.grid.MX(2)/2+1),:,:)=0;
        w_for_wave = real(ifft2(fft_for_wave));
        if  model.wave.advection.sto
            [model, sigma, coef_modulation] = ...
                fct_sigma(model,fft_for_wave);
%             if model.advection.sto % w is transported by sigma_adv dBt
%                 [~, sigma_wave,coef_modulation_sigma_wave] = ...
%                     fct_sigma(model,fft_w_local);
%             else % w is not transported by sigma_adv dBt
%                 [model, sigma_wave,coef_modulation_sigma_wave] = ...
%                     fct_sigma(model,fft_w_local);
%             if ~ model.advection.sto % w is not transported by sigma_adv dBt
%                 sigma = sigma_wave;
%             end
            if model.wave.advection.sto
                siz_mod = size(coef_modulation);
                bool_mod = ( prod(siz_mod(1:3)) == 1 );
                if ~ bool_mod
                    error(['Wave dynamics with heterogeneous ' ...
                        'time-uncorrelated velocity is not coded yet']);
                end
            end
        else
            sigma = 0;
        end
    else
        sigma = 0;
        w_for_wave = w ;
%         w_local = w;
%         sigma_wave = sigma;
    end
    if ~ model.wave.advection.sto
        w_local = w_for_wave ;
        for p =1:2
            grad_w_local(:,:,:,:,p) = ...
                real(ifft2(fct_grad(model,fft_for_wave(:,:,p,:))));
        end
    end
    
%     % Definition of currents for the wave dynamics
%     if ~ model.wave.advection.sto
%         fft_w_local = fft_w;
%         if model.wave.advection.smooth_LS_current
%             % Filtering
%             fft_w_local((model.wave.grid.MX(1)/2+1): ...
%                 (end-1-model.wave.grid.MX(1)/2+1),:,:,:)=0;
%             fft_w_local(:,(model.wave.grid.MX(2)/2+1): ...
%                 (end-1-model.wave.grid.MX(2)/2+1),:,:)=0;
%             w_local = real(ifft2(fft_w_local));
%         else
%             w_local = w;
%         end
%         for p =1:2
%             grad_w_local(:,:,:,:,p) = ...
%                 real(ifft2(fct_grad(model,fft_w_local(:,:,p,:))));
%         end
%     end
    
    % Finer time integration
    time_local = time;
    next_time_local = time;
    next_time = time + model.advection.dt_adv;
    while next_time_local < next_time
%     while time_local < next_time
        % Wave group without current
        w_group_0_forward = 0;
%         w_group_0_forward = fct_w_group(model,wave_group_forward(:,3:4,:,:));
        
        % Local dynamical time step
%         warning('Check size w_group_0 for CFL computation')
        model.wave.advection.dt_adv = fct_CFL(...
            model.grid, model.wave.advection,...
            model.sigma, permute(w_group_0_forward,[1 3 2 4]));
%         model.wave.advection.dt_adv
        time_local = next_time_local;
        next_time_local = time_local + model.wave.advection.dt_adv;
        if next_time_local > next_time
            model.wave.advection.dt_adv = next_time - time_local;
            next_time_local = next_time;
        end
        
        if  model.wave.advection.sto
            % Time-uncorrelated oceanic current velocity            
%             [model,sigma,coef_modulation,sigma_s,sigma_n_s] = ...
%                 fct_sigma(model,fft_w);
            sigma_dBt_on_sq_dt_local = ...
                fct_sample_sigma_dBt_on_sq_dt(model,...
                sigma,coef_modulation);
            % Full oceanic current velocity   
%             w_local = w ...
            w_local = w_for_wave ...
                + sigma_dBt_on_sq_dt_local / ...
                sqrt(model.wave.advection.dt_adv);
            for p =1:2
                grad_w_local(:,:,:,:,p) = ...
                    real(ifft2(fct_grad(model,fft2(w_local(:,:,p,:)))));
            end
%             grad_w_local = permute(grad_w_local,[1 2 5 4 3]);
            % Mx x My x 2(grad) x (N_pcl) x 2(velocity)
        end
        
        if model.wave.advection.sto
            siz_mod = size(coef_modulation);
            bool_mod = ( prod(siz_mod(1:3)) == 1 );
%             if ~(numel(coef_modulation) == 1)
            if ~ bool_mod
                error(['Wave dynamics with heterogeneous ' ...
                    'time-uncorrelated velocity is not coded yet']);
                heterogeneous_correction = bsxfun(@times, ...
                    coef_modulation, heterogeneous_correction) ;
                grad_w_local = grad_w_local + heterogeneous_correction;
            end
        end
                
        % Method of characteristics
        
        % Forward advection (to get k(x_0,t) on a regular grid)
        wave_group_forward = RK4_advection_lagrangienne(...
            model.wave.advection, model.grid, ...
            wave_group_forward, nan, ...
            w_local, grad_w_local, w_group_0_forward, X0);
    end
    
    % Check if there are Nan values
    if any(isnan(wave_group_forward(:)))
        error(['There are NaNs ' ...
            '(possibly because some advected points left the domain)']);
    end
    
    if model.advection.plot_charactericts
        % Save trajectory
%         ray_traject = cat(3,ray_traject,X_wave_forward);
%         ray_K = cat(3,ray_K,k_wave_forward);
%         ray_A = cat(3,ray_A,Ampli_wave_forward);
        ray_wave_group = cat(3,ray_wave_group,wave_group_forward);
    end
    
    
    %% Ray integration backward
%     % Finer time integration
%     time_local = time;
%     next_time_local = time;
%     next_time = time + model.advection.dt_adv;
%     
%     while next_time_local < next_time
% %     while time_local < next_time
%         % Wave group without current
%         w_group_0 = fct_w_group(model,wave_group(:,3:4,:,:));
%         
%         % Local dynamical time step
% %         warning('Check size w_group_0 for CFL computation')
%         model.wave.advection.dt_adv = fct_CFL(...
%             model.grid, model.wave.advection,...
%             model.sigma, permute(w_group_0,[1 3 2 4]));
% %         model.wave.advection.dt_adv
%         time_local = next_time_local;
%         next_time_local = time_local + model.wave.advection.dt_adv;
%         if next_time_local > next_time
%             model.wave.advection.dt_adv = next_time - time_local;
%             next_time_local = next_time;
%         end
%         
%         if  model.wave.advection.sto
%             % Time-uncorrelated oceanic current velocity            
% %             [model,sigma,coef_modulation,sigma_s,sigma_n_s] = ...
% %                 fct_sigma(model,fft_w);
%             sigma_dBt_on_sq_dt_local = ...
%                 fct_sample_sigma_dBt_on_sq_dt(model,sigma,coef_modulation);
%             % Full oceanic current velocity   
%             w_local = w ...
%                 + sigma_dBt_on_sq_dt_local / ...
%                 sqrt(model.wave.advection.dt_adv);
%             for p =1:2
%                 grad_w_local(:,:,:,:,p) = ...
%                     real(iff2(fct_grad(model,fft2(w_local(:,:,p,:)))));
%             end
% %             grad_w_local = permute(grad_w_local,[1 2 5 4 3]);
%             % Mx x My x 2(grad) x (N_pcl) x 2(velocity)
%         end
%         
%         % Method of characteristics
%         
%         if model.wave.advection.sto
%             if ~(numel(coef_modulation) == 1)
%                 error(['Wave dynamics with heterogeneous ' ...
%                     'time-uncorrelated velocity is not coded yet']);
%                 heterogeneous_correction = bsxfun(@times, ...
%                     coef_modulation, heterogeneous_correction) ;
%                 grad_w_local = grad_w_local + heterogeneous_correction;
%             end
%         end
%         
%         % Backward advection (to get k(x,t) on a regular grid)
%         wave_group = RK4_advection_lagrangienne(...
%             model.wave.advection, model.grid, ...
%             wave_group, ratio_ampli_ini, ...
%             - w_local, - grad_w_local, - w_group_0, X0);
%     end
%     
%     % Check if there are Nan values
%     if any(isnan(wave_group(:)))
%         error(['There are NaNs ' ...
%             '(possibly because some advected points left the domain)']);
%     end
    
    %% Adding time-correlated and time decorrelated velocity
    time = next_time;   
    
    if  model.advection.sto
        sigma_dBt_on_sq_dt = ...
            fct_sample_sigma_dBt_on_sq_dt(model,sigma_adv,coef_modulation);
        w = w + sigma_dBt_on_sq_dt/sqrt(model.advection.dt_adv);
    else
        sigma_dBt_on_sq_dt = 0;
    end
    % if isfield(model.advection, 'forcing') && model.advection.forcing.bool
    %         w(:,:,1) = w(:,:,1) + Vy;
    %     end
    
    %% Transport of tracer
    if ~ model.sigma.sto
        %if ~ model.sigma.sto
        % Runge-Kutta 4 scheme
        fft_b = RK4_fft_advection(model,fft_b, w);
    else
        % Euler scheme
        fft_b = fft_b ...
            + deriv_fft_advection( ...
            model, fft_b, w) ...
            * model.advection.dt_adv;
        clear model_temp
    end
    
    %% Discard particles which have blown up
    iii = isnan(fft_b) | isinf(abs(fft_b));
    if any(iii(:))
        iii=any(any(any(iii,3),2),1);
        if all(iii(:))
            if model.plots
                error('The simulation has blown up');
            else
                warning('One simulation has blown up');
                fprintf('Folder of the simulation : \n');
                fprintf([ model.folder.folder_simu ' \n']);
                return;
            end
        end
        nb_dead_pcl = sum(iii);
        warning([ num2str(nb_dead_pcl) ' particle(s) on ' num2str(N_ech) ...
            ' have(s) blown up and' ...
            ' are(is) resampled uniformly on' ...
            ' the set of the others particles']);
        N_ech_temp = N_ech - nb_dead_pcl;
        fft_b(:,:,:,iii)=[];
        iii_sample = randi(N_ech_temp,nb_dead_pcl,1);
        for k=1:nb_dead_pcl
            fft_b(:,:,:,end+1) = fft_b(:,:,:,iii_sample(k));
        end
    end
    clear iii
    
    %% Veloicity covariance and Eulerian absolute diffusivity
    if model.advection.cov_and_abs_diff
        cov_w(t) = 1/prod(model.grid.MX) * w_ref_cov * w(:) ;
    else
        cov_w = nan ;
    end
    
    %% Plot and save
    Plot_and_save
end
toc