function wave_group = init_Lagrangian_pcl(model)
% Intialize localization of the wave groups
%

%% Initial wave position
type_ray = 'center';

switch type_ray
    case 'bottom_boundary'
        X=nan(model.wave.init.n_ray,2);
        Lx = model.grid.dX(1)*(model.grid.MX(1)-1);
%         temp=(linspace(0,,n_ray+1))';
%         X(:,1) = temp(1:(end-1));
        X(:,1) = (Lx/model.wave.init.n_ray)*(0.5:model.wave.init.n_ray);
        X(:,2) = model.wave.init.X0;
    case 'center'
        Lx = model.grid.dX(1)*(model.grid.MX(1)-1);
        Ly = model.grid.dX(2)*(model.grid.MX(2)-1);
%         temp=(linspace(0,,n_ray+1))';
%         X(:,1) = temp(1:(end-1));
        x = (Lx/model.wave.init.n_ray)*(0.5:model.wave.init.n_ray);
        y = (Ly/model.wave.init.n_ray)*(0.5:model.wave.init.n_ray);
        [x,y] = ndgrid(x,y);
        X(:,1) = x(:);
        X(:,2) = y(:);
        clear x y
end

X = repmat(X,[1 1 1 model.wave.advection.N_ech]);

% %% Initilization of group wave-vectors
% K0 = nan(size(X));
% K0(:,1,:,:) = 0 ;
% K0(:,2,:,:) = model.wave.init.K0 ;
% 
% %% Initilization of group amplitudes
% A0 = model.wave.init.A0 * ones(size(X(:,1,:,:)));

%% Gather in one N-D array
wave_group(:,1:2,:,:) = X;
% wave_group(:,3:4,:,:) = K0;
% wave_group(:,5,:,:) = A0;

