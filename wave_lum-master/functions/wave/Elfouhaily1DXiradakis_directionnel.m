function [SPEC_bi, SPEC , Kp, ustar] = ...
    Elfouhaily1DXiradakis_directionnel(U,fetch,direction,grid)
%U wind speed in m/sec
%fetch in km
%kappa wavenumber in 1/m

direction_selection = true;

% grid
kappa = grid.k.k;
% m_MX = mean(grid.MX);
% phi = (0:(m_MX-1))*(2*pi/m_MX) ;
% phi_on_grid = atan2(grid.k.ky, grid.k.kx);
phi = atan2(grid.k.ky, grid.k.kx) - direction;

if direction_selection
    window = ( phi <= pi/2 ) & ( phi >= - pi/2 ) ;
    window = window / mean(window(:)) ;
else
    window = ones(grid.MX) ;
end

% Parameter of the bi_directionnal spectrum
% g=9.81;
km = 363;
a0=log(2)/4;
ap=4;
cm=0.23;

%JONSWAP 1-D Spectrum
% general costants
g = 9.81; %gravity acceleration in m/sec^2
Cm = 0.23; %capillary peak celerity in m/sec
Km = 370; %capillary peak wavenumber in 1/m
Xo = 22000; % fetch constant


X1 = fetch*1000; % dimensional fetch in m
% parameters
% low frequency curvature spectrum Bl calculation
% special costants
Ko = g/U.^2 ;
X = (g*X1)/U^2; % non dimensional fetch
omega = 0.84 .*(tanh((X/Xo).^0.4)).^(-0.75);
Kp = Ko.* omega.^2; % wavenumber of gravity wave peak
Ap = 0.006* sqrt(omega); % equilibrium range parameter
Cp = sqrt (g/Kp); % phase speed at spectral peak
c = sqrt (g./kappa); % phase speed
Lpm = exp((-5/4).*((Kp./kappa).^2)); % PM spectral shape
s = 0.08*(1+4.*omega^-3);
GAMA = exp(-((sqrt(kappa./Kp)-1).^2)/(2*s^2));

if omega <1
    gama = 1.7;
else gama = 1.7+log(omega);
end
Jp = gama.^GAMA; % JONSWAP peak enhancement factor
Z = exp( (- omega/sqrt(10)).*(sqrt(kappa./Kp)-1));
Bl = (Ap/2)*(Cp./c).*Lpm.* Jp.*Z; %low frequency curvature spectrum
%high frequency curvature spectrum Bh calculation
Zo = 3.7*(10^-5)*(U^2/g)*(U/Cp)^0.9; % roughness length
ustar = 0.42*U/log(10/Zo); %friction velocity
if ustar<Cm
    Am = 0.01*(1+log(ustar/Cm));
else Am= 0.01*(1+3*log(ustar/Cm));
end
Bh = (Am/2)*(Cm./c).*exp(-0.25*((kappa./Km)-1).^2).*Lpm; %high frequency curvature spectrum

SPEC = (Bh+Bl)./kappa.^3; % 1-D spectrum

%% Bi-directionnal spectrum computation

am=0.13*ustar/cm;

cp=sqrt (g/Kp);
c=sqrt( (g./kappa).*(1+(kappa./km).^2) );

% Spreading function
SF=1/(2*pi).*(1+cos(2.*phi).*(tanh(a0+ap.*(c./cp).^2.5+ am*(cm./c).^2.5) ));

% Bidirectionnal spectrum
SPEC_bi = window .* SPEC .* SF ./ kappa ;



