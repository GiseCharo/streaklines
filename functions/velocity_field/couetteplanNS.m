function dx=couetteplanNS(t,CIx,CIy)
%%couette plan non stationary computated over a grid of points
h=1;
U0=1;
nu=1;
omega=2*nu*((2*pi/h)^2);

% if nargin < 3
%     if length(CIx(:))~=2
%         error('Wrong size')
%     end
%     CIy = CIx(2);
%     CIx(2)=[];
% end

dx(:,:,1)=(U0/sinh(2*pi))*(-sin((2*pi/h)*CIy).*cosh((2*pi/h)*CIy).*sin(omega*t)...
    +cos((2*pi/h)*CIy).*sinh((2*pi/h)*CIy).*cos(omega*t));
dx(:,:,2)=0;
% dxx=(U0/sinh(2*pi))*(-sin((2*pi/h)*CIy).*cosh((2*pi/h)*CIy).*sin(omega*t)...
%     +cos((2*pi/h)*CIy).*sinh((2*pi/h)*CIy).*cos(omega*t));
% dxy=0;
% % if size(CI,1)==1
% %    dx=reshape(dx,[1 3 2]); 
% % end
% if nargin < 3
%     dx = [dxx ; dxy];
% else
%     dx(:,:,1)=dxx;
%     dx(:,:,2)=dxy;    
% end
 
end