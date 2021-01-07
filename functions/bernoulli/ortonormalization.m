function fn=ortonormalization(f)
normf=sqrt(f(:,1).*f(:,1)+f(:,2).*f(:,2));
fn=zeros(length(f),2);
for j=1:length(f)
fn(j,:)=f(j,:)./normf(j);
end
end