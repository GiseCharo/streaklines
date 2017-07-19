function mediaf=esperanza(f)

n=size(f,1);
if n>1
    m=mean(f);
    mediaf(:,1)=m(1,1)*ones(n,1);
    mediaf(:,2)=m(1,2)*ones(n,1);
else
    mediaf(1,1)=f(1,1);
    mediaf(1,2)=f(1,2);
end
end