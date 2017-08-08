function streakl=streakline_2(x0,time0,t,velocity)
% computes the streakline
%

streakl=zeros(length(time0),2);
%parfor i=1:length(time0)
for i=1:length(time0)
    streakl(i,:)=phi(x0,time0(i),t,velocity);
end
end