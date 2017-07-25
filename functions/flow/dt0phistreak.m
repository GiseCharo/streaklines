function Dt0phi=dt0phistreak(t0,x0,deltat0,t,velocity)
%%computes the the tangent along the streakline
Dt0phi=zeros(length(t0),2);
for j=1:length(t0)
    Dt0phi(j,:)=dtophi(x0,t0(j),t,deltat0,velocity);
end
end

