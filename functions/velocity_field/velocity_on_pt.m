function dx=velocity_on_pt(f_vel,t,CIx)
% Take the velocity on a single point
%

dx=f_vel(t,CIx(1),CIx(2));
dx = permute(dx,[3 1 2]);
 
end