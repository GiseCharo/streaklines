clear;clc;
%%Compute the correlation between the lagrangian velocity and the
%%correction of the tangent of the flux
x0(1,1)=1;x0(1,2)=0.5;
t=10;deltat0=0.0005;
t0=0:deltat0:t;
tau=0;
vel=@DGyreNS;
%vel=@DGyreS;
%vel=@couetteplanNS;
%vel=@couetteplanS;
for k=1:length(t0)-1
    
    if t0(k)+tau<=t
        Vlag(k,:)=lagrangianvelocity(x0,t0(k),t,vel);
        Dtiphi(k,:)=dtiphi(x0,t0(k),deltat0,t-t0(k),vel);
    end
    
end

meanVlag=esperanza(Vlag);
meanDtiphi=esperanza(Dtiphi);
corrVdtiphi=correlacion(Vlag,Dtiphi,meanVlag,meanDtiphi);

