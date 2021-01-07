function c=matrixcoef(dphi)
%%matriz of the coeficients of the projection of the streakline along 
%%the tangent direction(s), s ort
c=zeros(length(dphi),4,2);
s=ortonormalization(dphi(:,1,:));
sperp=ortogonal(s);
for j=1:size(dphi,2)
c(:,j,1)=projection_streak(dphi(:,j,:),s);
end
for j=1:size(dphi,2)
c(:,j,2)=projection_streak(dphi(:,j,:),sperp);
end

end