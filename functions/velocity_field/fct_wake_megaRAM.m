function dx=fct_wake_megaRAM(time,CIx,CIy)
% Compute the velocity of wake flow at Reynolds 100
%

persistent dt x y N U_global

%% First load
if isempty(U_global)
    fprintf('Please wait, initialisation of the velocity field.\n')
    % 1st load
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
if t < 1 || t > 10*N
    error('Outside the time interval')
end
U_local = permute(U_global(:,:,t,:),[1 2 4 3]);


%% Interpolation
s=size(CIx);
dx(:,1) =interp2(x,y,U_local(:,:,1)',CIx(:),CIy(:));
dx(:,2) =interp2(x,y,U_local(:,:,2)',CIx(:),CIy(:));
dx = reshape(dx,[s 2]);

