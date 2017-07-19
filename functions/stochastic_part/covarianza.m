function cov=covarianza(f)
%%calcula la covarianza de la funcion f, size(f)=n*n3*n4*n5
s=size(f);
n=s(1);%n
n2=prod(s(2:end));%n2=n3*n4*n5;
f=reshape(f,[n,n2]);
% mf=mean(f,1);
% f=bsxfun(@minus,f,mf);
mf=mean(f(:));
f=f-mf;
cov=zeros(floor(n/2),1);
for k=1:floor(n/2)
    temp=f(1:end-k+1,:).*f(k:end,:);
    cov(k)=(1/(n-k))*(sum(temp(:)));
    clear temp;
end
% function cov=covarianza(f,g,mediaf,mediag)
% n=size(f,1);
% f=f-mediaf;
% g=g-mediag;
% cov=(1/n)*(f(:,1)'*g(:,1)+f(:,2)'*g(:,2));
% end
