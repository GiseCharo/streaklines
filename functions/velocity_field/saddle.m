function dx=saddle(t,CIx,CIy)
%%compute the velocity field of the velocity field with a saddle point in
%%((A/lambda)*cos(omega*t),0)
omega=pi/5;
lambda=1;
A=0.1;
forcing=A*cos(omega*t);
if CIx<=2 && CIy<=1 && CIx>=0 && CIy>=0
    dx(:,:,1) = lambda*(CIx-1)-forcing;
    dx(:,:,2) = -lambda*CIy;
else
    dx(:,:,1) =0;
    dx(:,:,2) =0;
    
end

end