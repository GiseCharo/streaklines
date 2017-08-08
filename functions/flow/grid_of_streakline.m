function streakl=grid_of_streakline(X0,Y0,t0,deltat0,t,velocity)
% Computes a grid of streaklines
%

time0=t0:deltat0:t;
streakl=zeros([size(X0) 2 length(time0)]);
parfor i=1:length(time0)
    %for i=1:length(time0)
    streakltemp = nan([size(X0) 2]);
    [X,Y]=phi_grid(X0,Y0,time0(i),t,velocity);
    streakltemp(:,:,1)=X;
    streakltemp(:,:,2)=Y;
    streakl(:,:,:,i)=streakltemp;
    %     streakl(:,:,1,i)=X;
    %     streakl(:,:,2,i)=Y;
end
end