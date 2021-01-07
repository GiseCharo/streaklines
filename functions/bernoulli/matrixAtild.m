function Atild=matrixAtild(dphi)
%%the coeficients of this matrix are obtained from the taylor approximation
%%of order 4 of the streakline at time t_i
%%size(Atilde)=n*m*k, n=number of points of the streakline, m=number of 
%%unknowns= gradient f .s,gradient f .(s ort), laplacian f, (s ortog)'*Hess(f)*(s ortog),
%%(s)'*Hess(f)*(s ortog),
%%k=order of the taylor approximation of the function streakline
c=matrixcoef(dphi);

Atild=zeros(length(dphi),5,4);
Atild(:,1,1)=c(:,1,1);
Atild(:,1,2)=c(:,2,1)/2;
Atild(:,1,3)=c(:,3,1)/6;
Atild(:,1,4)=c(:,4,1)/24;
Atild(:,2,1)=c(:,1,2);
Atild(:,2,2)=c(:,2,2)/2;
Atild(:,2,3)=c(:,3,2)/6;
Atild(:,2,4)=c(:,4,2)/24;
Atild(:,3,2)=(c(:,1,1).^2)/2;
Atild(:,3,3)=(c(:,1,1).*c(:,2,1))/2;
Atild(:,3,4)=(c(:,2,1).^2)/8 +(c(:,1,1).*c(:,3,1))/6;
Atild(:,4,2)=(c(:,1,2).^2)/2- (c(:,1,1).^2)/2;
Atild(:,4,3)=(c(:,1,2).*c(:,2,2)-c(:,1,1).*c(:,2,1))/2;
Atild(:,4,4)=(c(:,2,2).^2)/8 -((c(:,2,1)).^2)/8+(c(:,1,2).*c(:,3,2)-...
    c(:,1,1).*c(:,3,1))/6;
Atild(:,5,2)=c(:,1,1).*c(:,1,2);
Atild(:,5,3)=(c(:,1,1).*c(:,2,2) +c(:,2,1).*c(:,1,2))/2;
Atild(:,5,4)=(c(:,2,1).*c(:,2,2))/4 +(c(:,1,1).*c(:,3,2)+c(:,1,2).*c(:,3,1))/6;
end