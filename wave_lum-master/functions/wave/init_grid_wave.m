function grid = init_grid_wave (grid)
% Create a grid in the spatial and Fourier space
%
% From struct "model", expecting:
%   - grid.MX (grid size);
%   - grid.dX (spatial sampling step);

grid.dX = grid.LX./grid.MX;
grid.x = grid.dX(1)*(0:(grid.MX(1)-1));
grid.y = grid.dX(2)*(0:(grid.MX(2)-1));
[grid.X,grid.Y]=ndgrid(grid.x,grid.y);

% check grid size is even
if any( mod(grid.MX,2)~=0)
    error('the number of grid points by axis need to be even');
end
PX = grid.MX/2;
ZM = PX + 1; %index of the single high-freq mode to be zero'ed out.

%% "Unstable" Fourier grid
% for homogeneous diffusion Laplacian(b), hyper-viscosity Laplacian^p(b),
% SQG relationship...
nx = [ 0:(PX(1)-1) 0 (1-PX(1)):-1]; %NB: the central single high-freq is zero'ed-out
ny = [ 0:(PX(2)-1) 0 (1-PX(2)):-1];
kx = (1./grid.MX(1)) .* nx;
ky = (1./grid.MX(2)) .* ny;
% the 2D grid
[kx,ky] = ndgrid(kx,ky);
kx = (2.*pi/grid.dX(1)) .* kx; %as wavenumbers
ky = (2.*pi/grid.dX(2)) .* ky;
k2 = kx.^2+ky.^2;
k2(ZM(1),:) = 0.; %de-alias the single high freq
k2(:,ZM(2)) = 0.;
k = sqrt(k2); %the modulus

%% Save
% the "unstable" grid
grid.k.ZM = ZM; %indices of the single mode to force to zero
grid.k.kx = kx;
grid.k.ky = ky;
grid.k.k2 = k2;
grid.k.k = k;

end
