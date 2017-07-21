function dx=fct_wake(time,CIx,CIy)
% Compute the velocity of wake flow at Reynolds 100
%

persistent dt x y U big_T_pers

%% First load
if isempty(dt)
    load([ pwd '/data/wakeflow/file_DNS100_inc3d_2017_07_17_1'],'dt','dX')
end
if isempty(big_T_pers)
    big_T_pers = - inf;
end

%% Chooose the right file
big_T = floor(time/(1000*dt))+1;
if big_T ~= big_T_pers
%     big_T
%     fprintf('Please wait, a new velocity file is loading.\n')
    if big_T < 1 || big_T > 20
        error('Outside the time interval')
    else
        load([ pwd '/data/wakeflow/file_DNS100_inc3d_2017_07_17_' ...
            num2str(big_T)],'U')
    end
    big_T_pers = big_T;
end

%% Extract velocity at the time of interest
t_local = floor(time/dt-1000*(big_T-1) )+1;
U_local = permute(U(:,:,t_local,:),[1 2 4 3]);

%% Grid
if isempty(x)
    [Mx,My,~]=size(U_local);
    x=(0:(Mx-1))*dX(1);
    y=(0:(My-1))*dX(2);
    y=y-mean(y);
end

%% Interpolation
s=size(CIx);
dx(:,1) =interp2(x,y,U_local(:,:,1)',CIx(:),CIy(:));
dx(:,2) =interp2(x,y,U_local(:,:,2)',CIx(:),CIy(:));
dx = reshape(dx,[s 2]);

