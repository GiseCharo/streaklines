function [H1,H2]=prebernoullieulerian(CIx,CIy,t,deltat0,nu,dX,velocity)
%%compute1H1=(1/2)*dtnorm2(U)-nu*<U,lapU>;H2=Uxw-dt(U)+nu*lapU
%%pint=scalar product of <U,laplacianoU>
U=velocity(t,CIx,CIy);%Mx My d

dtnormav2=dertnorm2(CIx,CIy,t,deltat0,velocity);

w=vort_mat(U,dX);%[Mx My d n2]
J=[0 -1; 1 0];
dtv=dertvel(CIx,CIy,t,deltat0,velocity);
UU=permute(U,[3 4 1 2]); % n1 n2 Mx My
lap=laplacian_mat(UU,dX);
lapU=permute(lap,[3 4 1 2]);
vxw=zeros(size(w,1),size(w,2),2);
pint=zeros(size(w,1),size(w,2));
for i=1:size(CIx,1)
    for j=1:size(CIx,2)
        vel(1,1)=U(i,j,1);
        vel(2,1)=U(i,j,2);
        lpu(1,1)=lapU(i,j,1);
        lpu(2,1)=lapU(i,j,2);
        vxw(i,j,:)=-w(i,j)*J*vel;
        pint(i,j)=vel'*lpu;
    end
    clear vel;
    clear lpu;
end

H1=(1/2)*dtnormav2-nu*pint;
H2=vxw-dtv+nu*lapU;
end
