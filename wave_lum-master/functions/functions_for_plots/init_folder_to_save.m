%% Init the folder paths and names to save the results

if model.advection.HV.bool
    add_subgrid_deter = ['_HV' '_' fct_num2str(model.advection.HV.order/2)];
elseif model.advection.Lap_visco.bool
    add_subgrid_deter = '_Lap_visco';
else
    add_subgrid_deter = '_no_deter_subgrid';
    %add_subgrid_deter = [];
end
if model.sigma.sto & model.sigma.assoc_diff
    add_subgrid_deter = [add_subgrid_deter '_assoc_diff'];
end
% if ( model.advection.HV.bool || model.advection.Lap_visco.bool) && ...
%         model.advection.Smag.bool
if ( ( model.advection.HV.bool | model.advection.Lap_visco.bool) & ...
        model.advection.Smag.bool ) | ...
        (model.sigma.sto & model.sigma.Smag.bool )
    % if ( model.advection.HV.bool | model.advection.Lap_visco.bool) & ...
    %  model.advection.Smag.bool
    add_subgrid_deter = [add_subgrid_deter '_Smag'];
    %     add_subgrid_deter = [add_subgrid_deter '_kappamax_on_kappad_' ...
    %         fct_num2str(model.advection.Smag.kappamax_on_kappad) ...
    %         '_dealias_ratio_mask_LS_' ...
    %         fct_num2str(model.grid.dealias_ratio_mask_LS)];
    if model.sigma.sto & model.sigma.Smag.bool & ...
            model.sigma.Smag.epsi_without_noise
        add_subgrid_deter = [add_subgrid_deter '_epsi_without_noise'];
    end
elseif model.sigma.sto & model.sigma.hetero_modulation
    add_subgrid_deter = [add_subgrid_deter '_hetero_modulation'];
    if  isfield(model.sigma,'hetero_energy_flux_prefilter') && ...
            model.sigma.hetero_energy_flux_prefilter
        add_subgrid_deter = [add_subgrid_deter ...
            '_prefilter'];        
    end
    if  isfield(model.sigma,'hetero_energy_flux_postfilter') && ...
            model.sigma.hetero_energy_flux_postfilter
        add_subgrid_deter = [add_subgrid_deter ...
            '_postfilter'];
    end
elseif model.sigma.sto & model.sigma.hetero_modulation_V2
    add_subgrid_deter = [add_subgrid_deter '_hetero_modulation_V2'];
    if  isfield(model.sigma,'hetero_energy_flux_prefilter') && ...
            model.sigma.hetero_energy_flux_prefilter
        add_subgrid_deter = [add_subgrid_deter ...
            '_prefilter'];        
    end
    if  isfield(model.sigma,'hetero_energy_flux_postfilter') && ...
            model.sigma.hetero_energy_flux_postfilter
        add_subgrid_deter = [add_subgrid_deter ...
            '_postfilter'];
    end
elseif model.sigma.sto & model.sigma.hetero_energy_flux
    add_subgrid_deter = [add_subgrid_deter '_hetero_energy_flux'];
    if model.sigma.hetero_energy_flux_v2
        add_subgrid_deter = [add_subgrid_deter '_v2'];
    end
    if model.sigma.hetero_energy_flux_averaging_after
        add_subgrid_deter = [add_subgrid_deter '_1_3rd_before_norm'];
    end
    if isfield(model.sigma,'kappa_VLS_on_kappa_LS')
        add_subgrid_deter = [add_subgrid_deter ...
            '_kappa_LS_on_kappa_VLS_' ...
            num2str(1/model.sigma.kappa_VLS_on_kappa_LS)];
    end
    if isfield(model.sigma,'kappaLSforEspi_on_kappamin')
        add_subgrid_deter = [add_subgrid_deter ...
            '_kappamin_on_kappaLSforEspi__' ...
            num2str(1/model.sigma.kappaLSforEspi_on_kappamin)];
    end
    if  isfield(model.sigma,'hetero_energy_flux_prefilter') && ...
            model.sigma.hetero_energy_flux_prefilter
        add_subgrid_deter = [add_subgrid_deter ...
            '_prefilter'];        
    end
    if  isfield(model.sigma,'hetero_energy_flux_postfilter') && ...
            model.sigma.hetero_energy_flux_postfilter
        add_subgrid_deter = [add_subgrid_deter ...
            '_postfilter'];
    end
elseif model.sigma.sto & model.sigma.hetero_modulation_Smag
    add_subgrid_deter = [add_subgrid_deter '_hetero_modulation_Smag'];
    if  isfield(model.sigma,'hetero_energy_flux_prefilter') && ...
            model.sigma.hetero_energy_flux_prefilter
        add_subgrid_deter = [add_subgrid_deter ...
            '_prefilter'];        
    end
    if  isfield(model.sigma,'hetero_energy_flux_postfilter') && ...
            model.sigma.hetero_energy_flux_postfilter
        add_subgrid_deter = [add_subgrid_deter ...
            '_postfilter'];
    end
end
if model.sigma.sto & model.sigma.no_noise
    add_subgrid_deter = [add_subgrid_deter '_no_noise'];
end

% if model.sigma.SelfSim_from_LS.bool
%     add_subgrid_deter = [add_subgrid_deter '_SelfSim_from_LS'];
% end

%%

if ~ model.sigma.sto % Deterministic case
    model.folder.folder_simu = [ 'images/usual_' model.dynamics ...
        add_subgrid_deter '/' model.type_data ];
elseif model.advection.sto &&  model.wave.advection.sto 
    % Stohastic wave % Full stochastic case
    %     model.folder.folder_simu = [ 'images/' model.dynamics ...
    %         '_MU' add_subgrid_deter '/' model.type_data ];
    model.folder.folder_simu = [ 'images/' model.dynamics ...
        '_MU' add_subgrid_deter '/' ...
        'type_spectrum_sigma_' model.sigma.type_spectrum '/' ...
        model.type_data ];
elseif model.wave.advection.sto % Stohastic wave
    %     model.folder.folder_simu = [ 'images/' model.dynamics ...
    %         '_MU' add_subgrid_deter '/' model.type_data ];
    model.folder.folder_simu = [ 'images/' model.dynamics ...
        '_Wave_LUM_' add_subgrid_deter '/' ...
        'type_spectrum_sigma_' model.sigma.type_spectrum '/' ...
        model.type_data ];
end
if model.wave.advection.smooth_LS_current
    model.folder.folder_simu = [ model.folder.folder_simu '_at_' ...
    num2str(model.wave.grid.MX(1)) 'x' num2str(model.wave.grid.MX(2)) ];
end
if model.advection.forcing.bool
    model.folder.folder_simu = [ model.folder.folder_simu ...
        '_forced_turb_' model.advection.forcing.forcing_type];
else
    model.folder.folder_simu = [ model.folder.folder_simu ...
        '_free_turb' ];
end
model.folder.folder_simu = [ model.folder.folder_simu ...
    '/' num2str(model.grid.MX(1)) 'x' num2str(model.grid.MX(2)) ];
if ( ( model.advection.HV.bool | model.advection.Lap_visco.bool) & ...
        model.advection.Smag.bool)
    subgrid_details = ['kappamax_on_kappad_' ...
        fct_num2str(model.advection.Smag.kappamax_on_kappad) ...
        '_dealias_ratio_mask_LS_' ...
        fct_num2str(model.advection.Smag.dealias_ratio_mask_LS)];
    if model.advection.Smag.spatial_scheme
        subgrid_details = [ subgrid_details '_spatial_scheme'];
    end
    model.folder.folder_simu = [ model.folder.folder_simu ...
        '/' subgrid_details ];
end
if model.sigma.sto
    if model.sigma.Smag.bool
        subgrid_details = ['kappamax_on_kappad_' ...
            fct_num2str(model.sigma.Smag.kappamax_on_kappad) ...
            '_dealias_ratio_mask_LS_' ...
            fct_num2str(model.advection.Smag.dealias_ratio_mask_LS)];
        if model.sigma.Smag.SS_vel_homo
            subgrid_details = [ subgrid_details '_SS_vel_homo'];
        elseif  model.sigma.proj_free_div
            subgrid_details = [ subgrid_details '_proj_free_div'];
        end
        if model.advection.Smag.spatial_scheme
            subgrid_details = [ subgrid_details '_spatial_scheme'];
        end
    elseif ( model.sigma.hetero_modulation ...
            |  model.sigma.hetero_modulation_V2 ...
            |  model.sigma.hetero_modulation_Smag ...
            | model.sigma.hetero_energy_flux )
        subgrid_details = ['dealias_ratio_mask_LS_' ...
            fct_num2str(model.advection.Smag.dealias_ratio_mask_LS)];
        if  model.sigma.proj_free_div
            subgrid_details = [ subgrid_details '_proj_free_div'];
        end
    end
    %     model.folder.folder_simu = [ model.folder.folder_simu ...
    %         '_kappamin_on_kappamax_' ....
    %         fct_num2str(model.sigma.kappamin_on_kappamax) ];
    if ~ ( exist('subgrid_details','var')==1)
        subgrid_details = [];
    end
    if ~ ( strcmp(model.sigma.type_spectrum,'EOF') || ...
            strcmp(model.sigma.type_spectrum,'Euler_EOF'))
        subgrid_details = [ subgrid_details ...
            '_kappamin_on_kappamax_' ....
            fct_num2str(model.sigma.kappamin_on_kappamax) ];
        if strcmp(model.sigma.type_spectrum,'Band_Pass_w_Slope')
            subgrid_details = [ subgrid_details ...
                '_on_kc_' ....
                fct_num2str(1/model.sigma.k_c) ];
        elseif strcmp(model.sigma.type_spectrum,'SelfSim_from_LS')
                if model.sigma.estim_k_LS
                    subgrid_details = [ subgrid_details ...
                     '_estim_k_LS'];
                end
                if model.sigma.time_smooth.bool
                    subgrid_details = [ subgrid_details ...
                        '_time_smooth_'... 
                        num2str(24*3600/model.sigma.time_smooth.tau)];
                end
        end
    else
        subgrid_details = [ subgrid_details ...
            '_nbDayLearn_' ...
            fct_num2str(model.sigma.nbDayLearn) ...
            '_Delta_T_on_Delta_t_' ...
            fct_num2str(model.sigma.Delta_T_on_Delta_t) ...;
            '_nb_EOF_' ...
            fct_num2str(model.sigma.nb_EOF)];
    end
    if model.sigma.N_ech > 1
        subgrid_details = [ subgrid_details ...
            '_N_ech_' ....
            fct_num2str(model.sigma.N_ech) ];
    end
    model.folder.folder_simu = [ model.folder.folder_simu ...
        '/' subgrid_details ];
end

% Create the folders
fct_create_folder_plots(model)

% Colormap
load('BuYlRd.mat');
model.folder.colormap = BuYlRd; clear BuYlRd

% Version of matlab
vers = version;
year = str2double(vers(end-5:end-2));
subvers = vers(end-1);
model.folder.colormap_freeze = ...
    (  year < 2014 || ( year == 2014 && strcmp(subvers,'a') ) );