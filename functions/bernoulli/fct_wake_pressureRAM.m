function pres=fct_wake_pressureRAM(time,CIx,CIy)
% Compute the pressure of wake flow at Reynolds 100
%

persistent dt x y N P_global

%% First load
if isempty(P_global)
    fprintf('Please wait, initialisation of the pressure field.\n')
    % 1st load
    load([ pwd '/data/wakeflow/file_DNS100_inc3d_2017_07_17_1'],...
        'dt','dX','P')
    % Grid
    [Mx,My,N]=size(P);
    x=(0:(Mx-1))*dX(1);
    y=(0:(My-1))*dX(2);
    y=y-mean(y);
    P_global = nan([Mx,My,10*N]);
    P_global(:,:,1:N) = P; clear P;
    for big_T = 2:10
        load([ pwd '/data/wakeflow/file_DNS100_inc3d_2017_07_17_' ...
            num2str(big_T)],'P')
        P_global(:,:, N*(big_T-1) + (1:N)) = P; clear P;
    end
    fprintf('Initialisation done.\n')
end

%% Extract velocity at the time of interest
t = floor(time/dt)+1;
if t < 1 || t > 10*N
    error('Outside the time interval')
end
P_local = permute(P_global(:,:,t),[1 2 4 3]);


%% Interpolation
s=size(CIx);
pres =interp2(x,y,P_local',CIx(:),CIy(:));
end
