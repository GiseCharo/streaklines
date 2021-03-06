function main(stochastic_simulation,type_data,resolution,forcing, ...
    sigma,Lap_visco,HV,Smag,advection_duration)
% Main function to Launch the code
%

%%%%%%%%%%%%%%%%%%%%
%%% Main
%%%%%%%%%%%%%%%%%%%%
if nargin == 0
    init;
end

% Choose to plot
% if nargin > 0
%     plots_bool = false;
% else
    plots_bool = true;
% end

%% Main parameters to choose

%% Waves

% Stochastic dispertion ratio
stochastic_wave = false

% Spatially-smooth the current for wave advection ?
smooth_LS_current = false;

% Resolution of the large_scale current
LS_resolution_wave = 128;
% LS_resolution_wave = 64;
    
% Periodic boundary conditions for wave group position
wave_periodic_boundary_conditions = true;

% Plot trajectories of wave groups
% (it may require a lot of RAM for many wave groups)
plot_charactericts = true;

% Initilization of group wave-vectors
wave.init.K0 = 2*pi/(3e1); % m^(-1) (wave length = 30 m)

% Initilization of group amplitudes
wave.init.A0 = 1e0; % m (wave heigth : 1 m )

% Initilization of group position
wave.init.X0 = 0; % m 
% wave.init.X0 = 1e4; % m 

% Number of pixels replicated on the border in order to deal with the
% boundaries, especially in the interpolation procedures
wave.advection.nbp=10;
% wave.advection.nbp=3;

% Number of rays
% wave.init.n_ray = 20;
wave.init.n_ray = 3;

%% Type of flow

% Type of dynamics
% dynamics = 'SQG';
dynamics = '2D';


bool_parfor = false;
bool_mat = true;
if bool_mat & bool_parfor
    error('No compatible');
end

if nargin == 0
    
    % Duration of the simulation (in seconds)
    advection_duration = 3600*24*30;
    % advection_duration = 3600*24*130;
    %advection_duration = 3600*24*1000;
%     advection_duration = 3600*24*16; % 20 days
%     advection_duration = 3600*24*102 % 20 days
    
%     bool_parfor = false
%     bool_mat = true
%     if bool_mat & bool_parfor
%         error('No compatible');
%     end
    
    % Type of initial condtions
    type_data ='disym_Vortices';
    % 'disym_Vortices' : 2 large dysymmetric  anticyclones and cyclones
    % 'Vortices' : 2 large anticyclones and 2 large cyclones
    %   (used in "Geophysical flow under location uncertainty", Resseguier V.,
    %    Memin E., Chapron B.)
    % 'Vortices2' : same as 'Vortices' but properly periodized (by Pierre Derian).
    % 'Perturbed_vortices' : Same flow with slight small-scale modifications
    %   (used in "Chaotic transitions and location uncertainty in geophysical
    %    fluid flows", Resseguier V., Memin E., Chapron B.)
    % 'Spectrum' : Gaussian random field with a spectrum slope deined by
    %   the variable slop_b_ini (default value  = -5/3)
    % 'Zero' : Field equal to zero everywhere
    % 'Constantin_case1'
    % 'Constantin_case2'
    
    % Resolution
%     resolution = 64
%         resolution = 128
    % resolution = 256
%     resolution = 512
%     resolution = 1024
    % resolution = 2048
    resolution = LS_resolution_wave
    
    % The number of grid point is resolution^2
    % It has to be an even integer
    
    % Forcing
    
    % Forcing or not
    forcing = true;
    % If yes, there is a forcing
    % F = ampli_forcing * odg_b * 1/T_caract * sin( 2 freq_f pi y/L_y)
    % % If yes, there is an additionnal velocity V = (0 Vy)
    % % with Vy = ampli_forcing * odg_b *  sin( 2 freq_f pi y/L_y)
end

% Type de forcing
% forcing_type = 'Kolmogorov';
forcing_type = 'Spring';

% Amplitude of the forcing
ampli_forcing = 30;
% ampli_forcing = 10;
% ampli_forcing = 1;

% Frequency of the forcing
freq_f = [3 2];
% freq_f = [0 1];

%%


if nargin == 0
    
    % Deterministic or random model
    stochastic_turb = false
    
    if stochastic_turb && ~ stochastic_wave
        warning('stochastic_turb is set to false');
        stochastic_turb = false
    end
    
    sigma.sto = stochastic_wave || stochastic_turb ;
    % Usual SQG model (stochastic_simulation=false)
    % or SQG_MU model (stochastic_simulation=true)
    
    if sigma.sto
        % Type of spectrum for sigma dBt
%         sigma.type_spectrum = 'Band_Pass_w_Slope'; % as in GAFD part II
        % sigma.type_spectrum = 'Low_Pass_w_Slope';
        % Spectrum cst for k<km ans slope for k>km
        % sigma.type_spectrum = 'Low_Pass_streamFct_w_Slope';
        % Matern covariance for the streamfunction
        % spectrum = cst. * k2 .* ( 1 + (k/km)^2 )^slope )
        % ~ k2 for k<km ans slope for k>km
        % sigma.type_spectrum = 'BB';
        % sigma.type_spectrum = 'Bidouille';
%         sigma.type_spectrum = 'EOF';
        %         sigma.type_spectrum = 'Euler_EOF';
        sigma.type_spectrum = 'SelfSim_from_LS'
        %  Sigma computed from self similarities from the large scales
        % sigma.type_spectrum = type_spectrum;
        
        % Homogeneous dissipation associated with the spectrum slope
        sigma.assoc_diff = false;
        
        % Smagorinsky-like control of dissipation
        sigma.Smag.bool = false;
        
        %     % Sigma computed from self similarities from the large scales
        %     sigma.SelfSim_from_LS.bool = true;
        
        %     if sigma.SelfSim_from_LS.bool
        %         % Sigma computed from a energy of absolute diffusivity spectrum
        %         % sigma.SelfSim_from_LS.spectrum = 'energy';
        %         sigma.SelfSim_from_LS.spectrum = 'abs_diff';
        %     end
        
        % if strcmp(sigma.type_spectrum,'SelfSim_from_LS')
        % Heterrogeenosu energy flux epsilon
        sigma.hetero_energy_flux = false;
        
        % if strcmp(sigma.type_spectrum,'SelfSim_from_LS')
        % Heterrogeenosu energy flux epsilon
        sigma.hetero_energy_flux_v2 = false;
        
%         if sigma.hetero_energy_flux_v2 
%             if ~ sigma.hetero_energy_flux
%                 error('sigma.hetero_energy_flux_v2 is always associated with sigma.hetero_energy_flux');
%             else
%                 sigma.hetero_energy_flux_averaging_after = true;
%             end
%         end
        
        % Modulation by local V L (estimated from the velocity and from
        % thegradient of the velocity)
        sigma.hetero_modulation = false;
        
        % Modulation by local V^2
        sigma.hetero_modulation_V2 = false;
        
        % Modulation Smag
        sigma.hetero_modulation_Smag = false;
        
        sigma.estim_k_LS = false;
        
        if strcmp(sigma.type_spectrum,'SelfSim_from_LS')
            %             % sigma.estim_k_LS = true;
            %             sigma.estim_k_LS = false;
            
            sigma.time_smooth.bool = false;
            %             sigma.time_smooth.bool = true;
            % sigma.time_smooth.tau = 24*3600 / 10;
            %             sigma.time_smooth.tau = 24*3600 / 2;
            sigma.time_smooth.tau = (64/resolution) * 24*3600 / 10 ;
        end
        
        if strcmp(sigma.type_spectrum,'EOF') ...
                || strcmp(sigma.type_spectrum,'Euler_EOF')
            % Ratio between the Shanon resolution and filtering frequency used to
            % filter the heterogenous diffusion coefficient
            Smag.dealias_ratio_mask_LS = 1/8;
            
            % Nb day used to learn to EOFs
            sigma.nbDayLearn= 50;
            % sigma.nbDayLearn= 500;
            
            % Time period sampling for the learning data set
            sigma.Delta_T_on_Delta_t = 8;
            % sigma.Delta_T_on_Delta_t = 4500;
            
            % Number of EOF (use all EOFs if set to inf)
            sigma.nb_EOF = 200; % ref
            % sigma.nb_EOF = 2;
            %             sigma.nb_EOF = 8000;
            % sigma.nb_EOF = inf;
        end
        %     %if strcmp(sigma.type_spectrum,'SelfSim_from_LS')
        %     if sigma.hetero_modulation & strcmp(sigma.type_spectrum,'SelfSim_from_LS')
        if sigma.hetero_modulation | sigma.hetero_energy_flux ...
                | sigma.hetero_modulation_V2 | sigma.hetero_modulation_Smag
            % Ratio between the Shanon resolution and filtering frequency used to
            % filter the heterogenous diffusion coefficient
            % Smag.dealias_ratio_mask_LS = 1/8;
%             Smag.dealias_ratio_mask_LS = 1/4;
%             Smag.dealias_ratio_mask_LS = 1/2;
            Smag.dealias_ratio_mask_LS = 1; % default value
            
            % Compute mudulation from filtered (kappa_min) fields
            sigma.hetero_energy_flux_prefilter = true;
            
            % Filter noise modulations (1/4*kappa_min) fields
            sigma.hetero_energy_flux_postfilter = true;
            
        end
        % end
        
        if sigma.hetero_energy_flux
            % Normalize before taking the power 1/3
            sigma.hetero_energy_flux_averaging_after = true;
            
            % To filter negative part of the energy flux
            
%             sigma.kappa_VLS_on_kappa_LS = sigma.kappamin_on_kappamax / 2
            sigma.kappa_VLS_on_kappa_LS = 1/8;
%             sigma.kappa_VLS_on_kappa_LS = 1/12;
%             sigma.kappa_VLS_on_kappa_LS = 1/16;
%             sigma.kappa_VLS_on_kappa_LS = 1/4;
%             sigma.kappa_VLS_on_kappa_LS = 1/4;
            
            % Maximum value considered 
            sigma.kappaLSforEspi_on_kappamin = 1;
        end
        
        if sigma.hetero_energy_flux_v2 
            if ~ sigma.hetero_energy_flux
                error('sigma.hetero_energy_flux_v2 is always associated with sigma.hetero_energy_flux');
%             else
%                 sigma.hetero_energy_flux_averaging_after = true;
            end
        end
        
        % Use a spatial derivation scheme for the herogeneous
        % disspation
        Smag.spatial_scheme = false;
        
        % Force sigma to be diveregence free
        sigma.proj_free_div = true;
        
        if ( (sigma.Smag.bool + sigma.hetero_modulation + ...
                sigma.hetero_energy_flux + sigma.hetero_modulation_V2 ...
                + sigma.hetero_modulation_Smag) > 1 ) ...
                || ( (sigma.Smag.bool + sigma.assoc_diff ) > 1 )
            error('These parametrizations cannot be combined');
        end
        
        if sigma.Smag.bool || sigma.assoc_diff
            % Rate between the smallest wave number of the spatially-unresolved
            % (not simulated) component of sigma dBt and the largest wave
            % number of the simulation
            sigma.kappaMinUnresolved_on_kappaShanon = 1;
            
            % Rate between the largest wave number of the spatially-unresolved
            % (not simulated) component of sigma dBt and the largest wave
            % number of the simulation
            sigma.kappaMaxUnresolved_on_kappaShanon = 8;
            
        end
        % For Smagorinsky-like diffusivity/viscosity or Hyper-viscosity,
        if sigma.Smag.bool
            % Smagorinsky energy budget (dissipation epsilon)
            % without taking into account the noise intake
            sigma.Smag.epsi_without_noise = false;
            
            % Ratio between the Shanon resolution and filtering frequency used to
            % filter the heterogenous diffusion coefficient
            % Smag.dealias_ratio_mask_LS = 1;
            Smag.dealias_ratio_mask_LS = 1/8;
            % Smag.dealias_ratio_mask_LS = 1/4;
            %Smag.dealias_ratio_mask_LS = 1/2;
            
            %         % Ratio between the Shanon resolution cut-off ( = pi / sqrt( dx*dy) )
            %         % and the targeted diffusion scale
            %         % %        sigma.Smag.kappamax_on_kappad = 2;
            %         % sigma.Smag.kappamax_on_kappad = 1;
            
            sigma.Smag.kappamax_on_kappad = 0.5; % (better(?))
            % sigma.Smag.kappamax_on_kappad = 1 / 4;
            %             sigma.Smag.kappamax_on_kappad = 1 / ...
            %                 sigma.kappaMaxUnresolved_on_kappaShanon;
            
            %         % Factor in front of the additional constant dissipation
            %         % Set to 0 for no additional constant dissipation
            %         sigma.Smag.weight_cst_dissip = 0;
            
            % Heterogeneity of the noise
            sigma.Smag.SS_vel_homo = false;
            
        end
        
        % Desactivate the noise
        sigma.no_noise = false;
        if sigma.no_noise
            warning('There is no noise here');
        end
    end
end

% Number of realizations in the ensemble
N_ech = 201
% N_ech = 200
% % ( N_ech=200 enables moments to converge when the parameter resolution is
% %   set to 128 )
% % ( N_ech is automatically set to 1 in deterministic simulations )
if ~sigma.sto
    N_ech = 1;
else
    N_ech
end


if nargin >0
    bool_parfor = false;
end

%% Deterministic subgrid tensor

if nargin == 0
    % Viscosity
    Lap_visco.bool = false;
    
    % % Smagorinsky-like viscosity
    % Smag.bool = false;
    % % HV.bool = false;
    
    % Hyper-viscosity
    HV.bool = true;
    
    if HV.bool
%         % HV.order=4;
%         HV.order=8;
        HV.order=2;
    end
    
    % Smagorinsky-like diffusivity/viscosity or Hyper-viscosity
    Smag.bool = false;
    
    % For Smagorinsky-like diffusivity/viscosity or Hyper-viscosity,
    if Smag.bool
        if Lap_visco.bool
            
            % Use a spatial derivation scheme for the herogeneous
            % disspation
            Smag.spatial_scheme = false;
            
            % Ratio between the Shanon resolution and filtering frequency used to
            % filter the heterogenous diffusion coefficient
            Smag.dealias_ratio_mask_LS = 1/8;
            %     dealias_ratio_mask_LS = 1/8;
            %     %dealias_ratio_mask_LS = 1/2;
            warning('Redondant argument that for heterogeneous small-scale velocity')
            
            % Ratio between the Shanon resolution cut-off ( = pi / sqrt( dx*dy) )
            % and the targeted diffusion scale
            %     %     Smag.kappamax_on_kappad = 2;
            %     % % Smag.kappamax_on_kappad = 1.5; % Better, still small oscillations or just pixels?
            %     % %    %  Smag.kappamax_on_kappad = 1.1; % Stable mais petit artefact
            %     Smag.kappamax_on_kappad = 1.1; % Stable mais petit artefact
            %     %  d'aliasing
            Smag.kappamax_on_kappad = 0.5;
            % Smag.kappamax_on_kappad = 1; % Stable mais petit artefact
            %  d'aliasing  % (better(?))
            
            % Factor in front of the additional constant dissipation
            % Set to 0 for no additional constant dissipation
            %    HV.weight_cst_dissip = 1/10;
            %     HV.weight_cst_dissip = 1/10; % no aliasing
            Smag.weight_cst_dissip = 0;
        elseif HV.bool
            % Ratio between the Shanon resolution and filtering frequency used to
            % filter the heterogenous diffusion coefficient
            Smag.dealias_ratio_mask_LS = 1;
            warning('Redondant argument that for heterogeneous small-scale velocity')
            
            % Ratio between the Shanon resolution cut-off ( = pi / sqrt( dx*dy) )
            % and the targeted diffusion scale
            Smag.kappamax_on_kappad = 1.1;% still small oscillations or just pixels?
            
            % Factor in front of the additional constant dissipation
            % Set to 0 for no additional constant dissipation
            % % %     HV.weight_cst_dissip = 1/10;% bit of (stable) aliasing
            %     % HV.weight_cst_dissip = 1/3; % still a bit of (stable) aliasing
            %     HV.weight_cst_dissip = 1/3;
            %     % HV.weight_cst_dissip = 0;
            
            Smag.weight_cst_dissip = 1/1;
            % % %     HV.weight_cst_dissip = 1/10;% bit of (stable) aliasing
        end
    else
        Smag.kappamax_on_kappad = 0;
    end
    %     % Smag.kappamax_on_kappad = 1.1;
    %     % % Smag.kappad_on_kappamax = 1/2;
    %     if Smag.bool
    %         warning('This value needed to be tuned?')
    %     end
end

%% Optional parameters

% Compute velocity covariance and absolute diffusivity
cov_and_abs_diff = false;

% Choose to plot one-point one-time moments each day
plot_moments = true;

% Choose to plot the dissipation by scale
plot_epsilon_k = false;
if sigma.sto & sigma.hetero_energy_flux
    plot_epsilon_k = true;
end

% Plot dissipations terms
plot_dissip = false;

% Begin simulation from a precomputed field?
use_save = false
% In this case, which day should be used as initialisation
% day_save = 15;
day_save = 100;
% day_save = 126;
if use_save
    day_save
end
%day_save = 150

dealias_method = 'exp';
% [WIP] Method for mandatory de-aliasing of non-linear terms in
% pseudospectral codes (advection and non-homogeneous stochastic diffusion)
% - 'lowpass': same as in SQGMU 1;
% - '2/3': the classical 2/3 rule (NB: produces Gibb's oscillations);
% - 'exp': high-order exponential filter (Constantin et al., J. Sci.
%   Comput. (2012)).

% Boundaries conditions
dirichlet = false;

% Variance tensor a_H
if sigma.sto
    switch sigma.type_spectrum
        case 'SelfSim_from_LS'
            sigma.k_c = 0;
        otherwise
            switch dynamics
                case 'SQG'
                    switch resolution
                        case 64
                            % sigma.k_c = 1/(3e2) /6 % 1/(1800 meters)
                            %  sigma.k_c = 1/(3e2) /10 % 1/(3000 meters)
                            sigma.k_c = 1/(3e2) /4 % 1/(1200 meters)
                            %                 % sigma.k_c = 1/(3e2) /3 % 1/(900 meters)
                        case 128
                            sigma.k_c = 1/(3e2) % 1/(300 meters)
                    end
                case '2D'
                    error(...
                        'The turbulence 2D is not stable under the action of noise');
                    %             k_c = 1/(eps); % 1/(100 meters)
                otherwise
                    error('Unknown type of dynamics');
            end
            % a_H is calculated latter in the code using
            % a_H = 2 * f_0 / k_c^2
            % where f_0 is the Corilis frequency
    end
else
    % If the simulation is deterministic, a_H = 0 and only one simulation
    % is performed
    sigma.k_c = inf; % And then a_H = 0
    N_ech=1;
    plot_moments = false;
end

% Spectrum slope of sigma dBt
switch dynamics
    case 'SQG'
        sigma.slope_sigma = - 5/3;
    case '2D'
        sigma.slope_sigma = - 3;
    otherwise
        error('Unknown type of dynamics');
end
sigma.slope_sigma_ref = sigma.slope_sigma;
% switch model.dynamics
%     case 'SQG'
%         sigma.slope_sigma_ref = -5/3;
%     case '2D'
%         sigma.slope_sigma_ref = -3;
%     otherwise
%         error('Unknown type of dynamics');
% end
if  sigma.sto & strcmp(sigma.type_spectrum,'BB')
    sigma.slope_sigma = 0;
    % elseif strcmp(sigma.type_spectrum,'SelfSim_from_LS')
    %     sigma.slope_sigma = nan;
end

if nargin == 0
    % Rate between the smallest and the largest wave number of sigma dBt
    if sigma.sto
        switch sigma.type_spectrum
            case {'SelfSim_from_LS','EOF','Euler_EOF'}
                
                pre_estim_slope=1e-1;
                %%
                pre_5 = 5e-2;
                sigma.kappamin_on_kappamax = ...
                    (log(1-pre_5)/log(pre_estim_slope))^(2/HV.order);
                %         sigma.kappamin_on_kappamax = ...
                %             (log(1-pre_estim_slope)/log(pre_estim_slope))^(2/HV.order);
                %%
                %                 pre=1e-2;
                %                 sigma.kappamin_on_kappamax = ...
                %                     (log(1-pre)/log(pre_estim_slope))^(2/HV.order);
                %%
                sigma.kappamin_on_kappamax_estim_slope = ...
                    (log(1-pre_estim_slope)/log(pre_estim_slope))...
                    ^(2/HV.order);
                
                sigma.kappaLS_on_kappamax = 1/8;
                
            otherwise
                switch resolution
                    case  128
                        sigma.kappamin_on_kappamax = 1/2;
                    case 64
                        sigma.kappamin_on_kappamax = 1/3;
                    otherwise
                        error('unknown');
                end
        end
        
        %         if strcmp(sigma.type_spectrum , 'SelfSim_from_LS')
        %             if sigma.Smag.bool | ...
        %                     Lap_visco.bool | ( HV.bool & (HV.order<=4) )
        %                 sigma.kappamin_on_kappamax = 1/2;
        %                 % sigma.kappamin_on_kappamax = 1/4;
        %                 % sigma.kappamin_on_kappamax = 1/8;
        %             elseif ( HV.bool & (HV.order==8) )
        %                 switch resolution
        %                     case  128
        %                 % sigma.kappamin_on_kappamax = 1/2;
        %                 pre=1e-2;
        %                 pre_estim_slope=1e-1;
        %                 sigma.kappamin_on_kappamax = ...
        %                     (log(1-pre)/log(pre_estim_slope))^(2/HV.order)
        %                 sigma.kappamin_on_kappamax_estim_slope = ...
        %                     (log(1-pre_estim_slope)/log(pre_estim_slope))...
        %                     ^(2/HV.order)
        %                     case 64
        %                 pre=1e-2;
        %                 pre_estim_slope=1e-1;
        %                 sigma.kappamin_on_kappamax = ...
        %                     (log(1-pre)/log(pre_estim_slope))^(2/HV.order)
        %                 sigma.kappamin_on_kappamax_estim_slope = ...
        %                     (log(1-pre_estim_slope)/log(pre_estim_slope))...
        %                     ^(2/HV.order)
        % %                 sigma.kappamin_on_kappamax = 0.45;
        % %                 % sigma.kappamin_on_kappamax = 1/3;
        %                     otherwise
        %                         error('unknown');
        %                 end
        %             else
        %                 warning('kappamin_on_kappamax may be inapropriate');
        %                 sigma.kappamin_on_kappamax = 1/2;
        %                 % sigma.kappamin_on_kappamax = 1/4;
        %                 % sigma.kappamin_on_kappamax = 1/8;
        %             end
        %
        %             sigma.kappaLS_on_kappamax = 1/8;
        %         else
        %             switch resolution
        %                 case  128
        %                     sigma.kappamin_on_kappamax = 1/2;
        %                 case 64
        %                     sigma.kappamin_on_kappamax = 1/3;
        %                 otherwise
        %                     error('unknown');
        %             end
        %
        % %             %kappamin_on_kappamax = 1/32;
        % %             sigma.kappamin_on_kappamax = 1/2;
        % %             % sigma.kappamin_on_kappamax = 1/128;
        % %             %         sigma.slope_sigma = - 5;
        % %             % warning('THIS PARAMETER NEEDS TO BE CHANGED -- TEST');
        %
        %             sigma.kappaLS_on_kappamax = 1/8;
        %         end
        
        % Rate between the largest wave number of sigma dBt and the largest wave
        % number of the simulation
        sigma.kappamax_on_kappaShanon = 1;
    end
end

% Spectrum slope of the initial condition (if type_data = 'Spectrum' )
switch dynamics
    case 'SQG'
%         slope_b_ini = 5/3;
        slope_b_ini = - 5/3;
    case '2D'
        slope_b_ini = - 3;
    otherwise
        error('Unknown type of dynamics');
end

% Physical parameters
model = fct_physical_param(dynamics);

% Gather parameters in the structure model
% model.wave.sigma = sigma;
% if stochastic_turb
% %     model.sigma.sto = sigma.sto;
%     model.sigma = sigma;
% end
model.sigma = sigma;
model.advection.sto = stochastic_turb;
if stochastic_turb
    model.advection.N_ech=N_ech;
else
    model.advection.N_ech=1;    
end
model.wave = wave;
model.wave.advection.sto = stochastic_wave;
if stochastic_wave
    model.wave.advection.N_ech=N_ech;
else
    model.wave.advection.N_ech=1;    
end
model.sigma.N_ech = max([ ...
    model.advection.N_ech model.wave.advection.N_ech]);
model.wave.advection.periodic_boundary_conditions = ...
    wave_periodic_boundary_conditions;
model.wave.grid.MX = LS_resolution_wave * [1 1];
model.wave.advection.smooth_LS_current = smooth_LS_current;
model.wave.advection.HV.bool = 0;
model.wave.advection.Smag.bool = 0;
if sigma.sto
    eval(['model.sigma.fct_tr_a = @(m,k1,k2) fct_norm_tr_a_theo_' ...
        model.sigma.type_spectrum '(m,k1,k2);']);
    % eval(['model.sigma.fct_tr_a = @(m,k1,k2,alpha) fct_norm_tr_a_theo_' ...
    %     model.sigma.type_spectrum '(m,k1,k2,alpha);']);
    % model.sigma.slope_sigma = slope_sigma;
    % model.sigma.kappamin_on_kappamax = kappamin_on_kappamax;
end
if strcmp(type_data,'Spectrum')
    model.slope_b_ini = slope_b_ini;
end
model.dynamics=dynamics;
model.type_data=type_data;
% model.advection.N_ech=N_ech;
% model.sigma.k_c = k_c;
model.advection.advection_duration=advection_duration;
model.advection.plot_charactericts = plot_charactericts;
model.advection.plot_epsilon_k = plot_epsilon_k;
model.advection.plot_dissip = plot_dissip;
model.advection.plot_moments = plot_moments;
model.advection.forcing.bool = forcing;
model.advection.forcing.ampli_forcing = ampli_forcing;
model.advection.forcing.freq_f = freq_f;
model.advection.forcing.forcing_type = forcing_type;
model.advection.HV = HV;
model.advection.cov_and_abs_diff = cov_and_abs_diff;
model.advection.Lap_visco = Lap_visco;
model.advection.Smag = Smag;
model.advection.use_save = use_save;
model.advection.day_save = day_save;
model.grid.dealias_method = dealias_method; %de-aliasing method
%model.Smag.dealias_ratio_mask_LS = dealias_ratio_mask_LS;
model.plots = plots_bool;

if model.sigma.sto
    if ( strcmp(model.sigma.type_spectrum,'EOF') || ...
            strcmp(model.sigma.type_spectrum,'Euler_EOF') ) ...
            && ( model.sigma.nb_EOF > model.advection.N_ech )
        warning(['The number of EOF is larger than the ensemble size.' ...
            ' Some EOFs are hence removed.']);
        model.sigma.nb_EOF = model.advection.N_ech;
    end
end

%% Random generator
rng('default'); %  The default settings are the Mersenne Twister with seed 0.

%% Generating initial buoyancy
[fft_buoy,model] = fct_buoyancy_init(model,resolution);

%% Advection
[fft_buoy_final, model] = fct_fft_advection_sto_mat(model, fft_buoy);
    
%% Re-plot
% fct_re_plot(model, fft_buoy);

%% Post-process plots

% last_day_plot_qq = 120;
% %resolution_HR = 512
% % resolution_HR = 1024
% % resolution_HR = 16 * resolution
% resolution_HR = 4 * resolution
% % resolution_HR = 8 * resolution
% if ~ use_save
%     first_day = 0;
%     day_save = 0;
% else
%     first_day = 100;
% %     first_day = 101;
%     % first_day = day_save;
%     day_save = 101;
% end
% 
% if plots_bool
%     post_process_error_grid(model.sigma.sto,...
%         model.type_data,resolution,resolution_HR,...
%         model.advection.forcing.bool,model.sigma,...
%         model.advection.Lap_visco,model.advection.HV,...
%         model.advection.Smag,model.advection.N_ech,...
%         first_day,advection_duration,day_save);
% %         first_day,last_day_plot_qq);
% end
% 
% resolution_HR