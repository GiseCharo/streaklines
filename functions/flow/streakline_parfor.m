function streakl=streakline_parfor(x0,t0,deltat0,t,velocity)
%%%computes the streakline
time0=t0:deltat0:t;
streakl=zeros(length(time0),2);
parfor i=1:length(time0)
%for i=1:length(time0)
    streakl(i,:)=phi(x0,time0(i),t,velocity);
end
end