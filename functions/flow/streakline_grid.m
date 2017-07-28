function streakl=streakline_grid(x0,t0,deltat0,t,velocity)
% Computes the streakline with the function phi_grid
%

time0=t0:deltat0:t;
streakl=zeros(length(time0),2);
x01=x0(1);
x02=x0(2);
parfor i=1:length(time0)
    [X,Y]=phi_grid(x01,x02,time0(i),t,velocity);
    streakl(i,:)=[X Y];
end
end