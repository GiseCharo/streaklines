function [H1,H2]=prebernoullieulerian(CIx,CIy,t,deltat0,nu,dX,velocity)
%%compute1H1=(1/2)*dtnorm2(U)-nu*<U,lapU>;H2=Uxw-dt(U)+nu*lapU
%%pint=scalar product of <U,laplacianoU>
U=velocity(t,CIx,CIy);%Mx My d

dtnormav2=dertnorm2(CIx,CIy,t,deltat0,velocity);

w=vort_mat(U,dX);%[Mx My d n2]
dtv=dertvel(CIx,CIy,t,deltat0,velocity);
UU=permute(U,[3 4 1 2]); % n1 n2 Mx My
lap=laplacian_mat(UU,dX);
lapU=permute(lap,[3 4 1 2]);
vxw=zeros(size(w,1),size(w,2),2);
pint=zeros(size(w,1),size(w,2));
vxw(:,:,1)=-w.*-U(:,:,2);
vxw(:,:,2)=-w.*U(:,:,1);
pint=U(:,:,1).*lapU(:,:,1)+U(:,:,2).*lapU(:,:,1);

H1=(1/2)*dtnormav2-nu*pint;
H2=vxw-dtv+nu*lapU;
end
