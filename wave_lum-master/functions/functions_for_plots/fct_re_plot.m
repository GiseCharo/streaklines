function fct_re_plot(model,  fft_b)
% Advection of buoyancy using SQG or SQG_MU model
%

tic
%% Folder to save plots and files
init_folder_to_save

%% Grid

% % Spatial grid
% model.grid.x = model.grid.dX(1)*(0:model.grid.MX(1)-1);
% model.grid.y = model.grid.dX(2)*(0:model.grid.MX(2)-1);
%
% % Grid in Fourier space
% model = init_grid_k (model);
%
% % Ensemble size
% % N_ech=model.advection.N_ech;
%
% %% Preprocessing of the Lagrangian particles
%
% % For now, we do not consider masks since we assume periodic boudary
% % conditions
%
% % Reference grid
% X0{1} = model.grid.dX(1)* ...
%     ((-model.wave.advection.nbp):...
%       (model.grid.MX(1)+model.wave.advection.nbp-1));
% X0{2} = model.grid.dX(2)* ...
%     ((-model.wave.advection.nbp):...
%       (model.grid.MX(2)+model.wave.advection.nbp-1));
% % X0{1} = repmat( model.grid.dX(1)* ((-nbp):(MX(1)+nbp-1)) ...
% %     ,[1 1 1 model.advection.N_ech]);
% % X0{2} = repmat( model.grid.dX(2)* ((-nbp):(MX(2)+nbp-1)) ...
% %     ,[1 1 1 model.advection.N_ech]);
%
% wave_group_ini = init_Lagrangian_pcl(model);
% wave_group = wave_group_ini;
% wave_group_forward = wave_group_ini;
% ratio_ampli_ini = wave_group_ini(:,5,:,:) ...
%    .* ( wave_group_ini(:,3,:,:).^2 + wave_group_ini(:,4,:,:).^2 ).^(-1/4) ;
%
% if model.advection.plot_charactericts
%     ray_wave_group = wave_group;
% end
%
% %% Initialisation of the spatial fields
%
% % Create several identical realizations of the intial buoyancy
% fft_b = repmat(fft_b(:,:,1),[1 1 1 model.advection.N_ech]);
%
% % Remove aliasing
% fft_b(model.grid.k.ZM(1),:,:,:)=0;
% fft_b(:,model.grid.k.ZM(2),:,:)=0;
%
% % Initial large-scale velocity
% fft_w = SQG_large_UQ(model, fft_b);
% w=real(ifft2(fft_w));
%
% % figure;imagesc(real(ifft2(fft_b))');axis xy;axis equal;colorbar;
% % figure;imagesc(w(:,:,1)');axis xy;axis equal;colorbar;
%
% % % Create several identical realizations of the intial buoyancy
% % fft_b = repmat(fft_b(:,:,1),[1 1 1 model.advection.N_ech]);
%
% %% Choice of the variance tensor a
% model = fct_init_sigma(model,fft_w);
%
% %% Forcing
% model = fct_set_forcing(model);
%
% %% Hyperviscosity
% model = fct_set_HV(model,w);
%
% %% Choice of time step : CFL
% model.advection.dt_adv = fct_CFL(model.grid,model.advection,model.sigma,w);
% % model.advection.dt_adv = fct_CFL(model,w);
%
% %% Statistics of the time-uncorrelated velocity
% [model,sigma,coef_modulation,sigma_s,sigma_n_s] = ...
%     fct_sigma(model,fft_w);

%% Loop on time

% Used for the first plot
tt_last = -inf;

% Number of time step
% N_t = ceil(model.advection.advection_duration/model.advection.dt_adv);

first_wave_plot = true;
fft_wave_noise = nan;

%% Print informations

% if model.plots
%     % Printing some information
%     fprintf(['The initial condition is ' model.type_data ' \n'])
%     fprintf(['1/k_c is equal to ' num2str(1/model.sigma.k_c) ' m \n'])
%     fprintf(['Time step : ' num2str(model.advection.dt_adv) ' seconds \n']);
%     fprintf(['Time of advection : ' num2str(...
%         model.advection.advection_duration/3600/24) ' days \n']);
%     fprintf(['Ensemble size : ' num2str(model.sigma.N_ech) ' realizations \n']);
%     fprintf(['Resolution : ' num2str(model.grid.MX(1)) ' x ' ...
%         num2str(model.grid.MX(2)) ' \n']);
%     if model.advection.HV.bool
%         str_subgridtensor = 'Hyper-viscosity';
%     elseif model.advection.Lap_visco.bool
%         str_subgridtensor = 'Laplacian diffusion';
%     else
%         str_subgridtensor = 'None';
%     end
%     if ( model.advection.HV.bool | model.advection.Lap_visco.bool) & ...
%             model.advection.Smag.bool
%         str_subgridtensor = ['Heterogneous ' str_subgridtensor];
%     end
%     fprintf(['Deterministic subgrid tensor : ' str_subgridtensor ' \n']);
%     fprintf(['Model type : ' add_subgrid_deter ' \n']);
%     if model.sigma.sto | model.advection.Smag.bool
%         fprintf(['Details : ' subgrid_details ' \n']);
%     end
%     if model.sigma.sto
%         fprintf(['type spectrum sigma :' model.sigma.type_spectrum ' \n']);
%     end
% end
%
% if model.sigma.sto & ...
%         strcmp(model.sigma.type_spectrum,'SelfSim_from_LS') & ...
%         model.sigma.time_smooth.bool
%     sigma_s = 0;
%     % sigma_s = sigma;
% end

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

bool_save = false;

while time < model.advection.advection_duration
    
    %     ratio_plot_by_day = 10;
    
    %     day_num = (floor(time/24/3600));
    %     idx_subday_num = time - day_num*24*3600;
    %     idx_subday_num = (floor(ratio_plot_by_day*idx_subday_num/24/3600));
    
    
    %     if day_num + idx_subday_num/ratio_plot_by_day > day_last_plot
    
    %% Plot and save
    Plot_and_save
    
    %%
    
    time = time + 24*3600/ratio_plot_by_day;
end
toc