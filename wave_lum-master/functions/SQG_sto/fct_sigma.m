function [model,sigma,coef_modulation,sigma_s,sigma_n_s] = ...
    fct_sigma(model,fft_w)
%% Statistics of the time-uncorrelated velocity
%

if model.sigma.sto
    if strcmp(model.sigma.type_spectrum,'SelfSim_from_LS')
        %             sigma = nan([model.grid.MX 2 N_ech]);
        %             tr_a = nan(1,N_ech);
        %             slope_sigma = nan(1,N_ech);
        %             offset_spectrum_a_sigma = nan(1,N_ech);
        %             km_LS = nan(1,N_ech);
        %%
        % for sampl=1:N_ech
        % parfor sampl=1:N_ech
        %%
        
        %             if model.sigma.time_smooth.bool
        %                 sigma_s = sigma;
        %             end
        
        % [sigma(:,:,:,sampl), ~, tr_a(sampl) ] ...
        [ sigma, ~, model.sigma.tr_a ,....
            model.sigma.slope_sigma,...
            model.sigma.offset_spectrum_a_sigma, ...
            model.sigma.km_LS ]...
            = fct_sigma_spectrum_abs_diff_mat(...
            model,fft_w,false);
        
        if model.sigma.time_smooth.bool
            sigma_n_s = sigma;
            d_sigma_s = 1/model.sigma.time_smooth.tau * ...
                ( - sigma_s + sigma_n_s ) ;
            sigma_s = sigma_s + d_sigma_s * model.advection.dt_adv;
            sigma = sigma_s;
        else
            sigma_s = nan;
            sigma_n_s = nan;
        end
        
        %             %                 % [sigma(:,:,:,sampl), ~, tr_a(sampl) ] ...
        %             %                 [ sigma(:,:,:,sampl), ~, tr_a(sampl) ,....
        %             %                     slope_sigma(sampl),...
        %             %                     offset_spectrum_a_sigma(sampl), ...
        %             %                     km_LS(sampl) ]...
        %             %                     = fct_sigma_spectrum_abs_diff(...
        %             %                     model,fft_w(:,:,:,sampl),false);
        model.sigma.a0 = model.sigma.tr_a/2;
        model.sigma.a0_on_dt = model.sigma.a0 / model.advection.dt_adv;
        % Diffusion coefficient
        model.advection.coef_diff = permute( 1/2 * model.sigma.a0 , ...
            [1 3 4 2]);
        if model.sigma.assoc_diff | model.sigma.Smag.bool
            % warning('deal with slope when there are several realizations')
            model.sigma.a0_SS = ...
                1/2 * fct_norm_tr_a_theo_Band_Pass_w_Slope(...
                model, ...
                model.sigma.kappaMinUnresolved_on_kappaShanon ...
                *(pi/sqrt(prod(model.grid.dX))), ...
                model.sigma.kappaMaxUnresolved_on_kappaShanon ...
                *(pi/sqrt(prod(model.grid.dX))));
            model.sigma.a0_LS = ...
                model.sigma.a0 ;
            model.sigma.a0 = ...
                model.sigma.a0 ...
                + model.sigma.a0_SS;
            model.sigma.a0_on_dt = model.sigma.a0 ...
                / model.advection.dt_adv;
            % Diffusion coefficient
            model.advection.coef_diff = permute(...
                1/2 * model.sigma.a0 ,[1 3 4 2]);
            %end
            
            
            if model.sigma.Smag.bool
                iii_non_degenerate_a0_SS = (model.sigma.a0_SS > eps);
                %if model.sigma.a0_SS > eps
                if any(iii_non_degenerate_a0_SS)
                    if model.sigma.Smag.epsi_without_noise
                        sigma(:,:,:,iii_non_degenerate_a0_SS) = ...
                            bsxfun(@times,...
                            permute( ...
                            sqrt(2./(...
                            model.sigma.a0_LS(iii_non_degenerate_a0_SS) ...
                            + model.sigma.a0_SS(iii_non_degenerate_a0_SS))) ...
                            , [ 1 3 4 2]) , ...
                            sigma(:,:,:,(iii_non_degenerate_a0_SS)) );
                    else
                        sigma(:,:,:,iii_non_degenerate_a0_SS) = ...
                            bsxfun( @times, ...
                            permute( ...
                            sqrt(2./...
                            model.sigma.a0_SS(iii_non_degenerate_a0_SS)) ...
                            , [ 1 3 4 2]) , ...
                            sigma(:,:,:,iii_non_degenerate_a0_SS) );
                    end
                    %                     else
                    %                         sigma = [];
                end
                if any(~iii_non_degenerate_a0_SS)
                    if strcmp(model.sigma.type_spectrum,'SelfSim_from_LS')
                        % The absolute diffusivity diagnosed from the large-scale
                        % kinematic spectrum is too weak. It suggests that there
                        % are few small scales and no subgrid terms is needed.
                        % Moreover, setting subgris terms to zero prevent numerical
                        % errors.
                        sigma(:,:,:,~iii_non_degenerate_a0_SS) = ...
                            zeros(size(sigma(:,:,:,~iii_non_degenerate_a0_SS)));
                        model.advection.coef_diff(~iii_non_degenerate_a0_SS) = 0;
                    else
                        error('Unknow case');
                    end
                end
                %                     if model.sigma.a0_SS > eps
                %                         if model_sampl(sampl).sigma.Smag.epsi_without_noise
                %                             sigma(:,:,:,sampl) = ...
                %                                 sqrt(2/(model_sampl(sampl).sigma.a0_LS ...
                %                                 + model_sampl(sampl).sigma.a0_SS)) ...
                %                                 * sigma(:,:,:,sampl);
                %                         else
                %                             sigma(:,:,:,sampl) = ...
                %                                 sqrt(2/model_sampl(sampl).sigma.a0_SS) ...
                %                                 * sigma(:,:,:,sampl);
                %                         end
                %                     elseif strcmp(model_sampl(sampl).sigma.type_spectrum,'SelfSim_from_LS')
                %                         % The absolute diffusivity diagnosed from the large-scale
                %                         % kinematic spectrum is too weak. It suggests that there
                %                         % are few small scales and no subgrid terms is needed.
                %                         % Moreover, setting subgris terms to zero prevent numerical
                %                         % errors.
                %                         sigma(:,:,:,sampl) = ...
                %                             zeros(size(sigma(:,:,:,sampl)));
                %                         model_sampl(sampl).advection.coef_diff = 0;
                %                     else
                %                         error('Unknow case');
                %                     end
            end
            
        end
    elseif model.sigma.Smag.bool | model.sigma.assoc_diff
        model.sigma.a0 = model.sigma.a0_LS + model.sigma.a0_SS;
    else
        % Variance tensor
        %parfor sampl=1:N_ech
        model.sigma.a0 = 2 * model.physical_constant.f0 ...
            / model.sigma.k_c^2;
        %end
    end
    
    
    if model.sigma.hetero_energy_flux
        
        %parfor sampl=1:N_ech
        % for sampl=1:N_ech
        
        coef_modulation = ...
            fct_epsilon_k_onLine(model,fft_b,fft_w);
        %end
    elseif model.sigma.hetero_modulation | ...
            model.sigma.hetero_modulation_V2
        if ( isfield(model.sigma,'hetero_energy_flux_prefilter')  ...
                &&    model.sigma.hetero_energy_flux_prefilter ) ...
                || (isfield(model.sigma,'hetero_energy_flux_postfilter')  ...
                &&    model.sigma.hetero_energy_flux_postfilter)
            error('not coded yet')
        end
        coef_modulation = ...
            fct_coef_estim_AbsDiff_heterogeneous(model,fft_w);
    elseif model.sigma.hetero_modulation_Smag
        % Heterogeneous dissipation coefficient
        if isfield(model.sigma,'hetero_energy_flux_prefilter')  ...
                &&    model.sigma.hetero_energy_flux_prefilter
            % Pre-filtering
            fft_b_for_modulation = bsxfun(@times, ...
                model.grid.k_aa_sigma_hetero_energy_flux.mask, ...
                fft_b);
        else
            fft_b_for_modulation = fft_b;
        end
        coef_modulation = fct_coef_diff(model,fft_b_for_modulation);
        clear fft_b_for_modulation
        
        if isfield(model.sigma,'hetero_energy_flux_postfilter')  ...
                &&    model.sigma.hetero_energy_flux_postfilter
            % Post-filtering
            coef_modulation = real(ifft2(bsxfun(@times, ...
                model.grid.k_aa_sigma_hetero_energy_flux_post.mask, ...
                fft2(coef_modulation))));
        end
        m_coef_modulation = mean(mean(coef_modulation,1),2);
        coef_modulation = bsxfun( @times, ...
            1./m_coef_modulation, coef_modulation);
        clear m_coef_modulation
    elseif model.sigma.Smag.bool
        % Heterogeneous dissipation coefficient
        coef_modulation = fct_coef_diff(model,fft_b);
        % Coefficient coef_Smag to target a specific diffusive scale
        coef_modulation = model.sigma.Smag.coef_Smag * coef_modulation ;
        
        %             %     figure(12);fct_spectrum( model,fft2(sigma_dBt_dt));
        %             %   figure(13);fct_spectrum( model,fft2(coef_diff_aa));
        %             %     figure(14);fct_spectrum( model,fft2(sqrt(coef_diff_aa)));
        %             %    figure(15);imagesc(sqrt(coef_diff_aa)');axis xy;axis equal
        %             %%
        %             % for sampl=1:N_ech
        %parfor sampl=1:N_ech
        
        % Coefficient coef_Smag to target a specific diffusive scale
        iii_non_degenerate_a0_SS = (model.sigma.a0_SS > eps);
        %if model_sampl(sampl).sigma.a0_SS > eps
        model.advection.coef_diff = zeros([1 1 1 model.advection.N_ech]);
        if model.sigma.Smag.epsi_without_noise
            model.advection.coef_diff(iii_non_degenerate_a0_SS) = 1;
        else
            % Taking into account the noise in the energy budget
            %                 model.advection.coef_diff(:,:,:,iii_non_degenerate_a0_SS) = ...
            model.advection.coef_diff(iii_non_degenerate_a0_SS) = ...
                permute( ...
                (1 + model.sigma.a0_LS(iii_non_degenerate_a0_SS) ./ ...
                model.sigma.a0_SS(iii_non_degenerate_a0_SS)) , ...
                [1 4 3 2] );
        end
        if any(~iii_non_degenerate_a0_SS)
            if strcmp(model.sigma.type_spectrum,'SelfSim_from_LS')
                % The absolute diffusivity diagnosed from the large-scale
                % kinematic spectrum is too weak. It suggests that there
                % are few small scales and no subgrid terms is needed.
                % Moreover, setting subgris terms to zero prevent numerical
                % errors.
                model.advection.coef_diff(~iii_non_degenerate_a0_SS) = 0;
            else
                error('Unknow case');
            end
        end
        %             %                 if model_sampl(sampl).sigma.a0_SS > eps
        %             %                     if model_sampl(sampl).sigma.Smag.epsi_without_noise
        %             %                         model_sampl(sampl).advection.coef_diff = 1;
        %             %                     else
        %             %                         % Taking into account the noise in the energy budget
        %             %                         model_sampl(sampl).advection.coef_diff = ...
        %             %                             (1 + model_sampl(sampl).sigma.a0_LS / ...
        %             %                             model_sampl(sampl).sigma.a0_SS) ;
        %             %                     end
        %             %                 elseif strcmp(model_sampl(sampl).sigma.type_spectrum,'SelfSim_from_LS')
        %             %                     % The absolute diffusivity diagnosed from the large-scale
        %             %                     % kinematic spectrum is too weak. It suggests that there
        %             %                     % are few small scales and no subgrid terms is needed.
        %             %                     % Moreover, setting subgris terms to zero prevent numerical
        %             %                     % errors.
        %             %                     model_sampl(sampl).advection.coef_diff = 0;
        %             %                 else
        %             %                     error('Unknow case');
        %             %                 end
        %             model.advection.coef_diff = ...
        %                 bsxfun( @times, ...
        %                 model.advection.coef_diff,...
        %                 coef_modulation) ;
        %             %end
        
        if model.sigma.Smag.SS_vel_homo
            coef_modulation = mean(mean(coef_modulation,2),1);
            %coef_modulation = mean(coef_modulation(:));
        end
    else
        coef_modulation = 1;
    end
    %% Variance tensor
    model.advection.coef_diff = ...
        bsxfun( @times, ...
        model.advection.coef_diff,...
        coef_modulation) ;
    %%
    if model.sigma.assoc_diff
        
        % for sampl=1:N_ech
        %parfor sampl=1:N_ech
        model.sigma.a0 = ...
            model.sigma.a0_LS ...
            + model.sigma.a0_SS;
        model.sigma.a0_on_dt = ...
            model.sigma.a0 / ...
            model.advection.dt_adv;
        % Diffusion coefficient
        %                 model.advection.coef_diff = ...
        %                     coef_modulation * ...
        %                     1/2 * model.sigma.a0;
        model.advection.coef_diff = bsxfun(@times, ...
            permute( 1/2 * model.sigma.a0 , [ 1 3 4 2]) , ...
            coef_modulation ) ;
        %end
    end
    
    
    
    
    % Maximum of the variance tensor
    %         % coef_modulation_a0 = max(max(coef_modulation,[],2),[],1);
    %         a0_temp = nan([1 1 1 N_ech]);
    a0_temp = nan([model.sigma.N_ech 1]);
    if size(coef_modulation,4)==1
        coef_modulation = repmat(coef_modulation,[1 1 1 model.sigma.N_ech]);
    end
    
    % for sampl=1:N_ech
    %parfor sampl=1:N_ech
    
    model.sigma.a0 = bsxfun(@times, ...
        permute(model.sigma.a0 , [ 1 3 4 2]) , ...
        coef_modulation ) ;
    
    
    
    if model.sigma.Smag.bool
        iii_non_degenerate_a0_SS = (model.sigma.a0_SS > eps);
        %if model_sampl(sampl).sigma.a0_SS > eps
        if model.sigma.Smag.epsi_without_noise
            model.sigma.a0(iii_non_degenerate_a0_SS) = ...
                model.sigma.a0(iii_non_degenerate_a0_SS) ...
                ./ (model.sigma.a0_SS(iii_non_degenerate_a0_SS) ...
                + model.sigma.a0_LS(iii_non_degenerate_a0_SS));
        else
            % Taking into account the noise in the energy budget
            model.sigma.a0(iii_non_degenerate_a0_SS) = ...
                model.sigma.a0(iii_non_degenerate_a0_SS) ...
                ./ model.sigma.a0_SS(iii_non_degenerate_a0_SS);
        end
        if any(~iii_non_degenerate_a0_SS)
            if strcmp(model.sigma.type_spectrum,'SelfSim_from_LS')
                % The absolute diffusivity diagnosed from the large-scale
                % kinematic spectrum is too weak. It suggests that there
                % are few small scales and no subgrid terms is needed.
                % Moreover, setting subgris terms to zero prevent numerical
                % errors.
                model.sigma.a0(~iii_non_degenerate_a0_SS) = 0;
            else
                error('Unknow case');
            end
        end
    end
    
    a0_temp = max(max( model.sigma.a0 ,[],2), [],1);
    model.sigma.a0 = max(a0_temp(:)) ; clear a0_temp;
else
    sigma = 0;
    model.sigma.a0 = 0;
    coef_modulation=nan;
    sigma_s=nan;
    sigma_n_s=nan;
    % sigma_on_sq_dt = 0;
end

