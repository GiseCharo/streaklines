function dx=dGyreS(t,x)    
epsilon=0.1;
amplitude=0.1;

dx = zeros(2,1);    % a column vector
a = epsilon;
b = 1 - 2*epsilon;
forcing = a*x(1).^2 + b*x(1);

dx(1) = -pi*amplitude*sin(pi*forcing).*cos(pi*x(2));
dx(2) = pi*amplitude*cos(pi*forcing).*sin(pi*x(2)).*(2*a*x(1) + b);

end