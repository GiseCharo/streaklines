function dx=fct_wake_megaRAM_2blocks(time,CIx,CIy)
% Compute the velocity of wake flow at Reynolds 100
%

persistent dt x y N U_global idx_block Mx My
meth = 'cubic';
% meth = 'linear';

%% First load
if isempty(U_global)
    fprintf('Please wait, initialisation of the velocity field.\n')
    % 1st load
    idx_block = 1;
    load([ pwd '/data/wakeflow/file_DNS100_inc3d_2017_07_17_1'],...
        'dt','dX','U')
    % Grid
    [Mx,My,N,d]=size(U);
    x=(0:(Mx-1))*dX(1);
    y=(0:(My-1))*dX(2);
    y=y-mean(y);
    U_global = nan([Mx,My,10*N,d]);
    U_global(:,:,1:N,:) = U; clear U;
    for big_T = 2:10
        load([ pwd '/data/wakeflow/file_DNS100_inc3d_2017_07_17_' ...
            num2str(big_T)],'U')
        U_global(:,:, N*(big_T-1) + (1:N),:) = U; clear U;
    end
    fprintf('Initialisation done.\n')
end

%% Extract velocity at the time of interest
t = floor(time/dt)+1;
idx = ceil(t/(10*N));
if t < 1 || idx > 2
    error('Outside the time interval')
elseif idx ~= idx_block
    idx_block = idx;
    fprintf('Please wait, reinitialisation of the velocity field.\n')
    U_global = nan([Mx,My,10*N,2]);
    for big_T = (1:10)
        load([ pwd '/data/wakeflow/file_DNS100_inc3d_2017_07_17_' ...
            num2str(10*(idx-1)+big_T)],'U')
        U_global(:,:, N*(big_T-1) + (1:N),:) = U; clear U;
    end
    fprintf('Reinitialisation done.\n')
end
t = t - (idx-1)*10*N;
U_local = permute(U_global(:,:,t,:),[1 2 4 3]);

%% Interpolation

s=size(CIx);

iii_outisde = CIx(:) <= 0 | abs(CIy(:)) > 6;
% iii_outside_unknown = CIx > 20;
dx(iii_outisde,1) = 1;
dx(iii_outisde,2) = 0;

dx(~iii_outisde,1) =interp2(x,y,U_local(:,:,1)',...
    CIx(~iii_outisde),CIy(~iii_outisde),meth);
dx(~iii_outisde,2) =interp2(x,y,U_local(:,:,2)',...
    CIx(~iii_outisde),CIy(~iii_outisde),meth);

dx = reshape(dx,[s 2]);

% %% Interpolation
% s=size(CIx);
% dx(:,1) =interp2(x,y,U_local(:,:,1)',CIx(:),CIy(:));
% dx(:,2) =interp2(x,y,U_local(:,:,2)',CIx(:),CIy(:));
% dx = reshape(dx,[s 2]);

