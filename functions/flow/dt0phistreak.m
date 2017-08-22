function Dt0phi=dt0phistreak(t0,x0,t,velocity,deltat0)
%%computes the the tangent along the streakline
Dt0phi=zeros(length(t0),2);
if nargin < 5
    deltat0 = t0(2)-t0(1);
end
for j=1:length(t0)
    Dt0phi(j,:)=dtophi(x0,t0(j),t,deltat0,velocity);
end
end

