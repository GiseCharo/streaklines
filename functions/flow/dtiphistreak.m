function Dtiphi=dtiphistreak(t0,x0,deltat0,t,velocity)
%%computes the correction of the tangent along the streakline
Dtiphi=zeros(length(t0),2);

for j=1:length(t0)-1
    Dtiphi(j,:)=dtiphi(x0,t0(j),deltat0,t-t0(j),velocity);
end
if t0(end)==t
    Dtiphi(end,:)=0;
else
    Dtiphi(end,:)=dtiphi(x0,t0(end),deltat0,t-t0(end),velocity);
end
end