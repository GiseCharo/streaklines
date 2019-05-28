function fft_noise = main_1dalle2...
    (time_global,X_wave,k_wave,Ampli_wave,...
    nb_fig,first_plot,fft_noise,...
    time_duration,title_)
% Main function to Launch the code
%

%%%%%%%%%%%%%%%%%%%%
%%% Main
%%%%%%%%%%%%%%%%%%%%
if nargin == 0
    init;
end

if nargin < 8
    time_duration = 1;
end

%% Main parameters to choose

% Grid
% % % % % Lx = 1e1; Ly = Lx;
% % % % Lx = 2e3; Ly = Lx;
% % % % Lx = 2e4; Ly = Lx;
% % % Lx = 3e2; Ly = Lx;
% % Lx = 1.5e2; Ly = Lx;
% Lx0 = 1.5e2; 
% Lx0 = 1e3; 
Lx0 = 5e1; 
Lx = 7e1; Ly = Lx;

% Mx = 64; My = Mx;
Mx = 128; My = Mx;
model_local.grid.MX=[Mx My]; model_local.grid.LX=[Lx Ly];
% model_local.grid.dX = model_local.grid.LX./model_local.grid.MX;
% model_local.grid.x = model_local.grid.dX(1)*(0:(model_local.grid.MX(1)-1));
% model_local.grid.y = model_local.grid.dX(2)*(0:(model_local.grid.MX(2)-1));
% model_local.grid.x = model_local.grid.x - X_wave(1);
% model_local.grid.y = model_local.grid.y - X_wave(2);
% % [model_local.grid.X,model_local.grid.Y]=meshgrid(model_local.grid.x,model_local.grid.y);
% [model_local.grid.X,model_local.grid.Y]=ndgrid(model_local.grid.x,model_local.grid.y);
model_local.grid = init_grid_wave (model_local.grid);
% model_local = init_grid_k (model_local);

model_local.time.dt = 0.1;
% % model_local.time.dt = 0.25;
% % % model_local.time.dt = 1;
% model_local.time.duration = 1;
% % model_local.time.duration = 10;
% % model_local.time.dt = 1;
% % model_local.time.duration = 50;
model_local.time.duration = time_duration;

% Physical constants
% Gravity
model_local.physical_constant.g = 9.81;

model_local.N_ech = 1;
% model_local.N_ech = 100;

model_local.dynamics = nan;

% Spectrum parameters
model_local.grid.wave_stat.spectrum_slope.mean = -4;
% model_local.grid.wave_stat.spectrum_slope.slope = 4/Lx;
model_local.grid.wave_stat.std = 1;
% % model_local.grid.wave_stat.corr_length = Lx/2;
% model_local.grid.wave_stat.corr_length = Lx;
model_local.grid.wave_stat.corr_length = Lx0;
% % model_local.grid.wave_stat.corr_length = 2*pi*1e2;
model_local.grid.wave_stat.type_spectrum = 'Constant'; % two plane wave
% model_local.grid.wave_stat.type_spectrum = 'Gaussian'; % two plane wave
% model_local.grid.wave_stat.type_spectrum = 'double_swell'; % two plane wave
% % model_local.grid.wave_stat.type_spectrum = 'Elfouhaily';
% model_local.grid.wave_stat.U10 = 10; % m.s-1
% model_local.grid.wave_stat.fetch = 80; % km
model_local.grid.wave_stat.direction = 0; % rad
% model_local.grid.wave_stat.direction = pi/4; % rad

% % model_local.folder.folder_simu = [ pwd '\images'];
% model_local.folder.folder_simu = [ pwd '/images/' ...
%     '1_patch/' ...
%     'local_resolution_' num2str(model_local.grid.MX(1)) '_' ...
%         num2str(model_local.grid.MX(2)) '_' ...
%     'domain_' num2str(model_local.grid.LX(1)) '_'  num2str(model_local.grid.LX(2)) '/' ...
%     'local_spectrum_' model_local.grid.wave_stat.type_spectrum ...
%     ];
% % model_local.folder.folder_simu = [ pwd '\images'];
% % mkdir(model_local.folder.folder_simu)
% mkdir(model_local.folder.folder_simu)

%% Creation of the initial spatial filter
fft_spatial_filter = fct_spatial_filter_wave(...
    model_local,model_local.grid,false,k_wave);
%     model_local,model_local.grid,true,k_wave);
fft_spatial_filter = Ampli_wave * fft_spatial_filter; 

% fft_spatial_filter = exp( 1i*( ...
%         X_wave(1)*model_local.grid.k.kx + ...
%         X_wave(2)*model_local.grid.k.ky )) ...
%         .* fft_spatial_filter; 
% % model_local.grid.X = model_local.grid.X + X_wave(1);
% % model_local.grid.Y = model_local.grid.Y + X_wave(2);

% fft_spatial_filter = fct_spatial_filter_wave(model_local);

% filter = ifft2(fft_spatial_filter);
% sum(abs(imag(filter(:)))) / sum(abs(real(filter(:))))

%% Intial wave field

% Fourier transform of white noise
% fft_noise = fft2( randn( [ model_local.grid.MX 1 model_local.N_ech]));
% fft_noise = fft2( randn( [ model_local.grid.MX 1 model_local.N_ech])) ...
%     /sqrt(prod(model_local.grid.MX));
if first_plot
    fft_noise = sqrt(prod(model_local.grid.MX));
%     fft_noise = fft2( randn( [ model_local.grid.MX 1 model_local.N_ech]));
% %     fft_noise = fft2( randn( [ model_local.grid.MX 1 model_local.N_ech])) ...
% %       /sqrt(prod(model_local.grid.MX));
% %     fft_noise = ( randn( [ model_local.grid.MX 1 model_local.N_ech]));
end
% Multiplication by the Fourier transform of the filter
fft_wave_ini = bsxfun(@times,fft_spatial_filter,fft_noise);
clear noise

% %% Intial wave field
% wave_ini = Ampli_wave * exp( 1i*( ...
%         k_wave(1)*model_local.grid.X + ...
%         k_wave(2)*model_local.grid.Y ));
% fft_wave_ini = fft2(wave_ini);

%% Creation of the temporal filter
omega = sqrt( model_local.physical_constant.g * model_local.grid.k.k );
fft_temporal_filter_time1 = exp( 1i * omega );

%% Intial wave field
for time = time_global : model_local.time.dt : ...
        (time_global+model_local.time.duration)
    % Temporal filter
%     fft_temporal_filter = cos( omega * time );
      fft_temporal_filter = fft_temporal_filter_time1 .^ time ;
    
%     temporal_filter = ifft2(fft_temporal_filter);
%     sum(abs(imag(temporal_filter(:)))) / sum(abs(real(temporal_filter(:))))

    % Multiplication by the Fourier transform of the filter
    fft_wave = bsxfun(@times, fft_temporal_filter, fft_wave_ini);
    clear noise
    
%     wave = ifft2(fft_wave);
%     sum(abs(imag(wave(:)))) / sum(abs(real(wave(:))))
    
    % Homogeneous wave
    wave = real(ifft2(fft_wave));
    % Plot
    plot_wave(model_local,wave(:,:,:,1),time,model_local.grid,...
        nb_fig,title_);
    
    % Test
% % %     wave = real(ifft2(fft_wave));
% %     wave =ifft2(fft_wave);
% %     std_wave = sqrt(mean(abs(wave(:)).^2))
% %     wave = real(wave);
%     std_real_wave = sqrt(mean(abs(wave(:)).^2))
% %     keyboard
end
