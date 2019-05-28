function fft_filter = fct_spatial_filter_wave(model,grid,bool_plot,k_wave)
% Compute the spatial filter which gives the satial correlation of the
% initial wave field
%

% direction_selection = true;
direction_selection = false;

if nargin < 3
    bool_plot = true;
end

if nargin < 2
    grid = model.grid;
end

grid_ref = grid;
if nargin >= 4
    kx = grid.k.kx - k_wave(1); 
    ky = grid.k.ky - k_wave(2);
    k2 = kx.^2+ky.^2;
    PX = model.grid.MX/2;
    ZM = PX + 1; %index of the single high-freq mode to be zero'ed out.
    k2(ZM(1),:) = 0.; %de-alias the single high freq
    k2(:,ZM(2)) = 0.;
    k = sqrt(k2); %the modulus
    grid.k.kx = kx;
    grid.k.ky = ky;
    grid.k.k2 = k2;
    grid.k.k = k;
end

PX=grid.MX/2;
k = grid.k.k;
d_kappa = 2*pi/sqrt(prod(grid.MX.* grid.dX));

%% 1D Spectrum

% % Largest wave number
% k_inf = kappa(min(PX));

switch grid.wave_stat.type_spectrum
    % switch model.wind_sea.wave_stat.type_spectrum
    
    % switch model.wind_sea.wave_stat.type_cov_patch
    
    case 'Constant'
        dkx = grid.k.kx(2,1) - grid.k.kx(1,1);
        dky = grid.k.ky(1,2) - grid.k.ky(1,1); 
        nx = (grid.k.kx/dkx <= 1/2 ) & (grid.k.kx/dkx > - 1/2 );
        ny = (grid.k.ky/dky <= 1/2 ) & (grid.k.ky/dky > - 1/2 );
        spectrum = nx.*ny ;
        spectrum = spectrum * prod(grid.MX) ;
        
        
%         spectrum = zeros(grid.MX) ;
%         spectrum(1,1) = prod(grid.MX)^2 ;
% %         spectrum(1,1) = prod(grid.MX) ;
        normalisation = false;
        omnidirectional = false;
        
    case 'Dirac'
        % spectrum = ones(grid.MX) ;
        spectrum = ones(grid.MX) * prod(grid.MX) ;
        normalisation = false;
        omnidirectional = false;
        
        
    case 'Gaussian'
        lcorr = grid.wave_stat.corr_length;
%         %         lcorr = model.wind_sea.wave_stat.slob_corr_len;
%         spectrum = 2*pi*lcorr^2 * ...
%         spectrum = 2*pi*(lcorr/mean(grid.LX))^2 * ...
        spectrum = (lcorr)^2 / (2*pi) * ...
            exp( - lcorr^2/2 * grid.k.k .^2 ) ;
        % From contiinuous to discrete space
        spectrum = (2*pi/grid.dX(1)) * (2*pi/grid.dX(2)) * spectrum ;
%         spectrum = spectrum / 
% %             exp( - 2*pi*lcorr^2 * grid.k.k .^2 ) ;
% %         spectrum = spectrum / prod(grid.LX) ;
%         spectrum = spectrum / prod(grid.MX) ;
% % %         spectrum = spectrum / sum(spectrum(:) );
% % %         spectrum = spectrum / ( sum(spectrum(:)) * prod(grid.MX) );
        grid.wave_stat.std = 1 ;
        normalisation = false;
%         normalisation = true;
%         omnidirectional = true;
        omnidirectional = false;
        
        direction_selection = true;
        
    case 'omnidirectional_Elfouhaily'
        U10 = 10; % m.s-1
        fetch = 80; % km
        %         U10 = 10; % m.s-1
        %         fetch = 80; % km
        spectrum = Elfouhaily1DXiradakis(U10,fetch,k) ;
        
        % Remove spatial mean
        spectrum(1,1,:,:)=0;
        
        normalisation = false;
        omnidirectional = true;
        
    case 'Elfouhaily'
        U10 = grid.wave_stat.U10; 
        fetch = grid.wave_stat.fetch; 
        direction = grid.wave_stat.direction;
%         U10 = 10; % m.s-1
%         fetch = 80; % km
%         direction = pi/4; % rad
        spectrum = ...
            Elfouhaily1DXiradakis_directionnel(U10,fetch,direction,grid);
        
%         figure();imagesc(spectrum'); axis xy; axis equal;colorbar;
%         keyboard;
        
        % Remove spatial mean
        spectrum(1,1,:,:)=0;
        
        normalisation = false;
        omnidirectional = false;
        
    case 'Self_similar'
        
        % Smallest wave number
        k0 = 2*pi/grid.wave_stat.corr_length;
        % k0 = 2*pi/model.wind_sea.wave_stat.corr_length;
        
        % Large-scale slope
        large_scale_slope = +10;
        
        spectrum = ( (k/k0 ) .^large_scale_slope ).* ...
            ( 1 + (k/k0 ).^2 ) .^ ...
            ( (grid.wave_stat.wave_stat.spectrum_slope.mean - large_scale_slope)/2 );
        %             ( (model.wind_sea.wave_stat.spectrum_slope.mean - large_scale_slope)/2 );
        
        % Remove spatial mean
        spectrum(1,1,:,:)=0;
        
        normalisation = true;
        omnidirectional = true;
        
    case 'cos'
        spectrum = zeros(grid.MX);
        spectrum(4,2) = 1;
        spectrum(end-2,end) = 1;
        normalisation = true;
        omnidirectional = false;
        
    case 'double_cos'
        spectrum = zeros(grid.MX);
        spectrum(4,2) = 1;
        spectrum(end-2,end) = 1;
        spectrum(7,3) = 1;
        spectrum(end-5,end-1) = 1;
        normalisation = true;
        omnidirectional = false;
        
    case 'swell'
        spectrum = zeros(grid.MX);
        spectrum(4,2) = 1;
        grid.wave_stat.std = grid.wave_stat.std * sqrt(2) ;
        normalisation = true;
        omnidirectional = false;
        
    case 'double_swell'
        spectrum = zeros(grid.MX);
        spectrum(4,2) = 1;
        spectrum(7,3) = 1;
        grid.wave_stat.std = grid.wave_stat.std * sqrt(2) ;
        normalisation = true;
        omnidirectional = false;
        
end

% Antialiasing
spectrum(PX(1)+1,:,:,:)=0;
spectrum(:,PX(2)+1,:,:)=0;

% figure;imagesc(spectrum');axis xy;axis equal


% if strcmp(grid.wave_stat.type_spectrum,'Gaussian')
%     spectrum = spectrum / ( sum(spectrum(:)) * prod(grid.MX) );
% end    

if normalisation
%     % Normalization and forced a given amplitude

%     %     spectrum = (model.wind_sea.wave_stat.std ^2) * ...
%     spectrum = (grid.wave_stat.std ^2) * ...
%         spectrum / ( sum(spectrum(:)) / prod(grid.MX) );
%     %         spectrum / ( sum(spectrum(:)) * d_kappa^2 / prod(grid.MX) );
%     % prod(grid.MX) = prod(grid.MX)^2/prod(grid.MX)
%     % is for the discrete Parseval * energy BY POINT
%     % / variance of the fft(spatial white noise)
    
    spectrum = (grid.wave_stat.std ^2) * ...
        spectrum / ( sum(spectrum(:)) / prod(grid.MX)^2 );
    % prod(grid.MX) = prod(grid.MX)^2
    % is for the discrete Parseval * energy BY POINT
end
% figure;imagesc(spectrum');axis xy;axis equal;
if omnidirectional
    % To remove the 2 pi which appear when we integrate the spectrum over the
    % wave-vector angles and add the (2 pi)^2 which appear when we go from k to
    % 2*pi*k
    % And discretisation to go from continuous to discrete Fourier transform
    % f_sigma = f_sigma * (2*pi);
    % f_sigma = 1/prod(grid.dX) * f_sigma;
    spectrum = (2*pi/prod(grid.dX)) * spectrum;
    
    % from omnidirectional spectrum to Fourier tranform square modulus
    % Division by k in dimension 2
    spectrum = bsxfun(@times, grid_ref.k.k.^(-1) , spectrum );
    % else
    
    % Antialiasing
    spectrum(PX(1)+1,:,:,:)=0;
    spectrum(:,PX(2)+1,:,:)=0;
    
    % Remove spatial mean
    spectrum(1,1,:,:)=0;
end

%% Directionnality

direction = grid.wave_stat.direction;

phi = atan2(grid.k.ky, grid.k.kx) - direction;

if direction_selection
    window = ( phi <= pi/2 ) & ( phi >= - pi/2 ) ;
%     window = window / mean(window(:)) ;
else
    window = ones(grid.MX) ;
end
spectrum = spectrum .* window;

% %% Wave group properties
% if nargin >= 4
%     spectrum = Ampli_wave * spectrum; 
% end
% if nargin >= 6
%     spectrum = exp( i( X_wave(1)*grid.k.kx + X_wave(2)*grid.k.ky )) ...
%         .* spectrum; 
% end

%%

% Division by prod(grid.MX) because of the variance white-in-space
% noise
% Multiplication by prod(grid.MX) because the spectrum respresent
% an energy by unit of space and not the sum over the space
% f_sigma = 1/prod(grid.MX) * f_sigma;
% f_sigma = prod(grid.MX) * f_sigma;

% % Influence of discretisation
% f_sigma = f_sigma * d_kappa ;
% f_sigma = f_sigma / ( prod(MX.*dX) /(2*pi) ) ;

% % Influence of discretisation (from discrete to continuous spectrum)
% spectrum = spectrum * d_kappa^2 ;

% From square modulus to modulus
fft_filter = sqrt( spectrum );

% % From 1D function to 2D function
% f_sigma = interp1(kappa,[zeros(1,model.advection.N_ech); f_sigma],k);
% %f_sigma = interp1(kappa,[0; f_sigma],k);

% Cleaning
% f_sigma(k<=k0,:)=0;
% fft_filter(k>k_inf,:)=0;
% f_sigma=reshape(f_sigma,[grid.MX 1 model.advection.N_ech]);


if bool_plot
    id_part = 1;
    
    %% Wave number
    M_kappa=min(grid_ref.MX);
    P_kappa= M_kappa/2;
    kappa= d_kappa * ( 0:(P_kappa-1) ) ;
    
    %% Masks associated with the rings of iso wave number
    k = grid_ref.k.k;
    k=k(:);
    kappa_shift = [ kappa(2:end) kappa(end)+d_kappa ];
    idx = sparse( bsxfun(@le,kappa, k ) );
    idx = idx & sparse( bsxfun(@lt,k, kappa_shift ) );
    
    % Choose one realization
    spectrum_estim = fft_filter(:,:,:,id_part);
    
    % Compute the (discrete) bi-directional spectrum of sigma dBt
    spectrum_estim=abs(spectrum_estim).^2;
    
    % Compute the (discrete) omnidirectional spectrum of sigma dBt
    spectrum_estim = idx' * spectrum_estim(:);
    
    % Influence of the complex brownian variance
    spectrum_estim = prod(grid.MX)*spectrum_estim;
    
    % Division by prod(grid.MX) because of the Parseval theorem for
    % discrete Fourier transform
    % Division by prod(grid.MX) again in order to the integration
    % of the spectrum over the wave number yields the energy of the
    % height averaged (not just integrated) over the space
    spectrum_estim = 1/prod(grid.MX)^2 * spectrum_estim;
    
    %%
    
    % Division by the wave number step (=> continuous spectrum)
    spectrum_estim = spectrum_estim / d_kappa;
    
    % Wave std (for test)
    std_wave_estim = sqrt( d_kappa * sum(spectrum_estim(:)) ) ;
    
    
    %% Plot spectrum
    
    taille_police = 12;
    X0 = [0 0];
    
    figure1=figure(18);
    widthtemp = 12;
    heighttemp = 6;
    set(figure1,'Units','inches', ...
        'Position',[X0 widthtemp heighttemp], ...
        'PaperPositionMode','auto');
    loglog(kappa(2:end),spectrum_estim(2:end),'r.-')
    ax=axis;
    ax(4)=max(spectrum_estim);
    ax(4)=ax(4)*2;
    ax(3)=(kappa(2)/kappa(end))*min(max(spectrum_estim));
    ax(3) = min( [ax(3) min(spectrum_estim)]);
    ax(1:2)=kappa([2 end]);
    if ax(4)>0
        axis(ax)
    end
    ax = axis;
    hold off
    set(gca,'XGrid','on','XTickMode','manual');
    width = 9;
    height = 3;
    set(figure1,'Units','inches', ...
        'Position',[X0 width height], ...
        'PaperPositionMode','auto');
    set(gca,'YGrid','on')
    set(gca,...
        'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontSize',taille_police,...
        'FontName','Times')
    ylabel('$S$',...
        'FontUnits','points',...
        'interpreter','latex',...
        'FontSize',taille_police,...
        'FontName','Times')
    xlabel('$\kappa \bigl ( rad.m^{-1} \bigr )$',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontSize',taille_police,...
        'interpreter','latex',...
        'FontName','Times')
    title('Wave spectrum',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'interpreter','latex',...
        'FontSize',12,...
        'FontName','Times')
    
    %% Save plot
    drawnow
    folder_simu = model.folder.folder_simu;
%     mkdir(folder_simu)
    eval( ['print -depsc ' folder_simu '\wave_spectrum.eps']);
    
    std_wave_estim
end
end
