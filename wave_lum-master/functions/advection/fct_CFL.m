function dt = fct_CFL(grid,advection,sigma,w)
% Compute the CFL
%


% CFL of the diffusion (or CFL of the white noise advection)
dX2=(grid.dX /pi).^2;
if advection.sto
    % if sigma.sto
    bound1=2/sigma.a0*prod(dX2)/sum(dX2);
else
    bound1 = inf;
end

% CFL of the (large-scale) advection
dX=permute(grid.dX,[1 3 2]);
bound2=sum(bsxfun(@times,abs(w),pi./dX),3);
bound2=max(bound2(:));
bound2=1/bound2/4;

if advection.HV.bool
    % CFL of the hyperviscosity
    bound3=1/advection.HV.maxVal*(prod(dX2)/sum(dX2)) ^ ...
        (advection.HV.order/2);
else
    bound3 = inf;
end
clear dX dX2

if isfield(advection, 'forcing') && advection.forcing.bool
    bound4 = (1/advection.forcing.on_T)/100;
else
    bound4= inf;
end

% Minimum of the CFL
dt = min([bound1 bound2 bound3 bound4]);
clear bound1 bound2 bound3 bound4
if advection.sto
    dt=dt/2;
    % Further constraint on dt due to the use of a (simple) Euler scheme
    % for the SPDE
end
if advection.Smag.bool & advection.Smag.spatial_scheme
    dt=dt/2;
end
