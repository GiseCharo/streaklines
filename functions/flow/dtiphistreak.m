function Dtiphi=dtiphistreak(t0,x0,deltat0,t,velocity)
%%computes the correction of the tangent along the streakline
Dtiphi=zeros(length(t0),2);
for j=1:length(t0)-1
    Dtiphi(j,:)=dtiphi(x0,t0(j),deltat0,t-t0(j),velocity);
end
Dtiphi(end,:)=(phi(x0,t0(end),t0(end)+t-t0(end),velocity)-phi(x0,t0(end-1)...
    ,t0(end-1)+t-t0(end-1),velocity))/deltat0;
end

