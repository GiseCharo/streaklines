function [H1streak,H2streak]=interpstreak(stline,H1,H2,X,Y,t0)

%%interpolates the H1,H2 along the streakline

H1interp= griddedInterpolant(X,Y,H1(:,:),'spline');
H2xinterp= griddedInterpolant(X,Y,H2(:,:,1),'spline');
H2yinterp= griddedInterpolant(X,Y,H2(:,:,2),'spline');

H1streak=zeros(length(t0),1);
H2streak=zeros(length(t0),2);
H1streak=H1interp(stline(:,1),stline(:,2));
H2streak(:,1)=H2xinterp(stline(:,1),stline(:,2));
H2streak(:,2)=H2yinterp(stline(:,1),stline(:,2));
end