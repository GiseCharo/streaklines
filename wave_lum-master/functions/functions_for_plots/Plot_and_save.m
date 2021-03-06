
%% Plots and save
%     % tt = floor(t ); % Number of days
%     tt = floor(t *model.advection.dt_adv/ (3600*24)); % Number of days
%     tt = floor(time/ (3600*24)); % Number of days
%     if tt > tt_last
%         tt_last = tt;
%         day = num2str(floor(time/24/3600));
%         % day = num2str(floor(t*model.advection.dt_adv/24/3600));
%         % t_last_plot = t;

ratio_plot_by_day = 1;

day_num = (floor(time/24/3600));
idx_subday_num = time - day_num*24*3600;
idx_subday_num = (floor(ratio_plot_by_day*idx_subday_num/24/3600));

%%
%     warning('DEBUG')
%     % if (t_loop - t_last_plot)*dt >= 3600*24*1
%     day_num = (floor(time/24/3600));
%     day = num2str(day_num);
%     day
%     day_last_plot = day_num;
%
%     nb=2;
%     fft_w=fft_w(:,:,:,1:nb);
%     model.advection.N_ech =nb;
%     fct_sigma_spectrum_abs_diff_mat(model,fft_w,true,day);
%%


if day_num + idx_subday_num/ratio_plot_by_day > day_last_plot 
% if day_num > day_last_plot
    % if (t_loop - t_last_plot)*dt >= 3600*24*1
    day_num = (floor(time/24/3600));
    idx_subday_num = time - day_num*24*3600;
    idx_subday_num = (floor(ratio_plot_by_day*idx_subday_num/24/3600));

%     day_num = (floor(time/24/3600));
    day = num2str(day_num);
    idx_subday = num2str(idx_subday_num);
    day = [day '_' idx_subday];
    day
    day_last_plot = day_num + idx_subday_num/ratio_plot_by_day;
%     day_last_plot = day_num;


    %% Load files
    if ~ bool_save
        load( [model.folder.folder_simu '/files/' day '.mat']);
        fft_w = fft2(w);
    end
    
    if model.plots
        toc
        tic
        %             if ~model.sigma.sto
        %                 model_sampl = model;
        %             end
        
        fprintf([ num2str(time/(24*3600)) ' days of advection \n'])
        if model.sigma.sto
            a_0_LS = mean(sigma_dBt_on_sq_dt(:).^2);
            a_0_LS
        end
        
        id_part=1;
        %%
        if model.advection.Smag.bool || ...
                (model.sigma.sto && ( model.sigma.Smag.bool ...
                || model.sigma.hetero_modulation ...
                || model.sigma.hetero_modulation_V2 ...
                || model.sigma.hetero_energy_flux ...
                || model.sigma.hetero_modulation_Smag ) )
            % Coefficient coef_Smag to target a specific diffusive scale
            if model.advection.Smag.bool
                [coef_diff_aa,coef_diff] = fct_coef_diff(model,...
                    fft_b(:,:,:,id_part));
                coef_diff_aa = ...
                    model.advection.Smag.coef_Smag * coef_diff_aa ;
                coef_diff = ...
                    model.advection.Smag.coef_Smag * coef_diff ;
            elseif model.sigma.hetero_modulation_Smag
                if isfield(model.sigma,'hetero_energy_flux_prefilter')  ...
                        &&    model.sigma.hetero_energy_flux_prefilter
                    % Pre-filtering
                    fft_b_for_modulation = bsxfun(@times, ...
                        model.grid.k_aa_sigma_hetero_energy_flux.mask, ...
                        fft_b);
                else
                    fft_b_for_modulation = fft_b;
                end
                [coef_diff_aa,coef_diff] = fct_coef_diff(model,...
                    fft_b_for_modulation(:,:,:,id_part));
                
                
                if isfield(model.sigma,'hetero_energy_flux_postfilter')  ...
                        &&    model.sigma.hetero_energy_flux_postfilter
                    % Post-filtering
                    coef_modulation = real(ifft2(bsxfun(@times, ...
                        model.grid.k_aa_sigma_hetero_energy_flux_post.mask, ...
                        fft2(coef_modulation))));
                    coef_diff_aa = real(ifft2(bsxfun(@times, ...
                        model.grid.k_aa_sigma_hetero_energy_flux_post.mask, ...
                        fft2(coef_diff_aa))));
                end
                
                coef_diff_aa = coef_diff_aa / mean(coef_diff_aa(:));
                coef_diff = coef_diff / mean(coef_diff(:));
            elseif model.sigma.Smag.bool
                [coef_diff_aa,coef_diff] = fct_coef_diff(model,...
                    fft_b(:,:,:,id_part));
                if model.sigma.a0_SS(id_part) > eps
                    coef_diff = ...
                        (1 + model.sigma.a0_LS(id_part) ...
                        / model.sigma.a0_SS(id_part)) * ...
                        model.sigma.Smag.coef_Smag * coef_diff ;
                    coef_diff_aa = ...
                        (1 + model.sigma.a0_LS(id_part) ...
                        / model.sigma.a0_SS(id_part)) * ...
                        model.sigma.Smag.coef_Smag * coef_diff_aa ;
                elseif strcmp(model.sigma.type_spectrum,'SelfSim_from_LS')
                    % The absolute diffusivity diagnosed from the large-scale
                    % kinematic spectrum is too weak. It suggests that there
                    % are few small scales and no subgrid terms is needed.
                    % Moreover, setting subgris terms to zero prevent numerical
                    % errors.
                    coef_diff = 0;
                    coef_diff_aa = 0;
                else
                    error('Unknow case');
                end
                
            elseif model.sigma.hetero_modulation ...
                    | model.sigma.hetero_modulation_V2
                if isfield(model.sigma,'hetero_energy_flux_prefilter')  ...
                        &&    model.sigma.hetero_energy_flux_prefilter
                    % Pre-filtering
                    fft_b_for_modulation = bsxfun(@times, ...
                        model.grid.k_aa_sigma_hetero_energy_flux.mask, ...
                        fft_b);
                else
                    fft_b_for_modulation = fft_b;
                end
                
                fft_w_for_modulation = SQG_large_UQ(model, fft_b_for_modulation);
                %                     coef_diff_aa = ...
                %                         model.sigma.a0(id_part)/2  ;
                [coef_diff_aa,coef_diff] = ...
                    fct_coef_estim_AbsDiff_heterogeneous(...
                    model,fft_w_for_modulation(:,:,:,id_part));
                
                if isfield(model.sigma,'hetero_energy_flux_postfilter')  ...
                        &&    model.sigma.hetero_energy_flux_postfilter
                    % Post-filtering
                    coef_modulation = real(ifft2(bsxfun(@times, ...
                        model.grid.k_aa_sigma_hetero_energy_flux_post.mask, ...
                        fft2(coef_modulation))));
                    coef_diff_aa = real(ifft2(bsxfun(@times, ...
                        model.grid.k_aa_sigma_hetero_energy_flux_post.mask, ...
                        fft2(coef_diff_aa))));
                end
                
                coef_diff_aa = ...
                    model.sigma.a0(id_part)/2 * coef_diff_aa ;
                coef_diff = ...
                    model.sigma.a0(id_part)/2 * coef_diff ;
            elseif model.sigma.hetero_energy_flux
                %                     coef_diff_aa = ...
                %                         model.sigma.a0(id_part)/2  ;
                coef_diff_aa = ...
                    fct_epsilon_k_onLine(model,fft_b,fft_w);
                coef_diff_aa = coef_diff_aa(:,:,:,id_part);
                coef_diff_aa = ...
                    model.sigma.a0(id_part)/2 * coef_diff_aa ;
                coef_diff=nan(size(coef_diff_aa));
            end
            figure(9);
            subplot(1,2,1);
            imagesc(model.grid.x_ref,model.grid.y_ref,...
                real(ifft2( coef_diff))');axis xy;axis equal;colorbar;
            subplot(1,2,2);
            imagesc(model.grid.x_ref,model.grid.y_ref,...
                coef_diff_aa');axis xy;axis equal;colorbar;
            drawnow
            eval( ['print -depsc ' model.folder.folder_simu ...
                '/dissip_coef/' day '.eps']);
            %%
            figure(9);
            imagesc(model.grid.x_ref,model.grid.y_ref,...
                coef_diff_aa');axis xy;axis equal;colorbar;
            drawnow
            eval( ['print -depsc ' model.folder.folder_simu ...
                '/dissip_coef/' day '.eps']);
        end
        %    if model.advection.Smag.bool || model.sigma.Smag.bool ...
        %              || model.sigma.hetero_modulation
        %  % Coefficient coef_Smag to target a specific diffusive scale
        %         if model.sigma.Smag.bool
        %  [coef_diff_aa,coef_diff] = fct_coef_diff(model,fft_b);
        %   coef_diff_aa = model.sigma.Smag.coef_Smag * coef_diff_aa ;
        %          else
        %   coef_diff_aa = model.advection.Smag.coef_Smag * coef_diff_aa ;
        %          end
        %         figure(9);
        %        subplot(1,2,1);
        %         imagesc(model.grid.x_ref,model.grid.y_ref,...
        %         real(ifft2( coef_diff))');axis xy;axis equal;colorbar;
        %                 subplot(1,2,2);
        %                 imagesc(model.grid.x_ref,model.grid.y_ref,...
        %                     coef_diff_aa');axis xy;axis equal;colorbar;
        %                 drawnow
        %        eval( ['print -depsc ' model.folder.folder_simu ...
        %                     '/dissip_coef/' day '.eps']);
        %             end
        %%
        % Plots
        [spectrum,name_plot,int_epsilon] = ...
            fct_plot(model,fft_b,day,...
            ray_wave_group(:,1:2,:,:),1);
%             ray_wave_group(:,1:2,:,:),ray_wave_group(:,5,:,:));
        
        if model.advection.plot_dissip
            fct_plot_dissipation(model,fft_b,sigma,day);
        end
        
        if model.sigma.sto & ...
                strcmp(model.sigma.type_spectrum,'SelfSim_from_LS')
            fct_sigma_spectrum_abs_diff_mat(model,fft_w,true,day);
            %    sigma_on_sq_dt = (1/sqrt(model.advection.dt_adv)) ...
            %      * sigma; clear sigma
            %     model.sigma.a0 = a0;
            %    model.sigma.a0_on_dt = ...
            %        model.sigma.a0 / model.advection.dt_adv;
            %     % Diffusion coefficient
            %    model.advection.coef_diff = 1/2 * model.sigma.a0;
            %        warning('The CFL should be changed');
        end
        
        
        if model.sigma.sto & ...
                ( model.sigma.Smag.bool | model.sigma.assoc_diff )
            slope_sigma=model.sigma.slope_sigma(id_part)
            a0_LS=model.sigma.a0_LS(id_part)
            a0_SS=model.sigma.a0_SS(id_part)
        end
        
        if model.advection.cov_and_abs_diff
            abs_diff = sum(cov_w(t_ref_cov:end))*model.advection.dt_adv;
            figure(36)
            plot(model.advection.dt_adv*(0:(N_t-1))/3600/24,cov_w);
            hold on;
            plot(t_ref_cov*[1 1]*model.advection.dt_adv/3600/24,...
                max(abs(cov_w(~isnan(cov_w))))*[-1 1],'r');
            hold off
        end
        
        % Dissipation by scale
        if model.advection.plot_epsilon_k
            fct_plot_epsilon_k(model,fft_b(:,:,:,id_part),day);
            % fct_plot_epsilon_k(model,fft_b,int_epsilon,day);
        end
        dt = model.advection.dt_adv
        
    end
    
    %% Plot local wave groups
    first_wave_plot = true;
    dt_plot_wave = 0.2;
%     dt_plot_wave = 0.1;
    time_duration_plot_wave = 2 ;
    iidx_plot_wave = -1;
%     for time_plot_wave = time:dt_plot_wave:(time+time_duration_plot_wave)
%         iidx_plot_wave = iidx_plot_wave + 1;
%         for p=1:3
% %             wave_group_local = mean(wave_group_forward(p,:,:,:),4);
%             wave_group_ampli_local_med = ...
%                 median(wave_group_forward(p,5,:,:),4);
%             iii_wave_group_ampli_local_med = ...
%                 find(squeeze(wave_group_ampli_local_med == ... 
%                   wave_group_forward(p,5,:,:) ));
%             iii_wave_group_ampli_local_med = ...
%                     iii_wave_group_ampli_local_med(1);
%             wave_group_local = wave_group_forward(p,:,:, ...
%                 iii_wave_group_ampli_local_med);
%             
% %             fft_wave_noise = main_1dalle(time_plot_wave,...
%             fft_wave_noise = main_1dalle2(time_plot_wave,...
%                 wave_group_local(:,1:2,:),...
%                 wave_group_local(:,3:4,:),...
%                 wave_group_local(:,5,:),...
%                 p,first_wave_plot,fft_wave_noise, ...
%                 0, num2str(day_last_plot));
% %             fft_wave_noise = main_1dalle(time_plot_wave,...
% %                 wave_group_forward(p,1:2,:,id_part),...
% %                 wave_group_forward(p,3:4,:,id_part),...
% %                 wave_group_forward(p,5,:,id_part),...
% %                 p,first_wave_plot,fft_wave_noise, ...
% %                 0, num2str(day_last_plot));
%             first_wave_plot = false;
%             %         wave_group_forward(p,1:2,:,id_part)
%             %         wave_group_forward(p,3:4,:,id_part)
%         end
%         figure12=figure(12);drawnow;
%         eval( ['print -dpng ' model.folder.folder_simu '/local_wave/'  ...
%             num2str(day) '_' num2str(iidx_plot_wave) '.png']);
%         figure4=figure(4);drawnow;
%         eval( ['print -dpng ' model.folder.folder_simu '/local_wave_spectrum/'  ...
%             num2str(day) '_' num2str(iidx_plot_wave) '.png']);
%     end
    
    % eval( ['print -dpng ' folder_simu '/Spectrum_' name_save '_' num2str(time) '.png']);
    
    %% Save files
%     if bool_save
%         save( [model.folder.folder_simu '/files/' ...
%             num2str(day) '_' num2str(iidx_plot_wave) ...
%             '.mat'], ...
%             'model','time','fft_b','w','sigma_dBt_on_sq_dt', ...
%             'sigma','wave_group', 'wave_group_forward');
%         if model.advection.plot_charactericts
%             save( [model.folder.folder_simu '/files/'  ...
%             num2str(day) '_' num2str(iidx_plot_wave) ...
%             '.mat'], ...
%                 'ray_wave_group','-append');
%         end
%     end
    if bool_save
        save( [model.folder.folder_simu '/files/' day '.mat'], ...
            'model','time','fft_b','w','sigma_dBt_on_sq_dt', ...
            'sigma','wave_group', 'wave_group_forward');
        if model.advection.plot_charactericts
            save( [model.folder.folder_simu '/files/' day '.mat'], ...
                'ray_wave_group','-append');
        end
    end
    
    %             'sigma_on_sq_dt');
    %         %             'model','t','fft_b','w','sigma_dBt_on_sq_dt', ...
    %         %             'sigma_on_sq_dt');
    %         %         %             'sigma_on_sq_dt','cov_w','abs_diff');
end

%     % Dissipation by scale
%     if (t == t_last_plot + 1) &  model.advection.plot_epsilon_k
%         fct_plot_epsilon_k(model,fft_b,day);
%         % fct_plot_epsilon_k(model,fft_b,int_epsilon,day);
%     end
