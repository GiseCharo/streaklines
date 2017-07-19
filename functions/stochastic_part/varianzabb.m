function sigma=varianzabb(f,deltat0)
%%f es el incremento de la funcion
%%f(t0+deltat0)-f(t0)
n=size(f,1);
sigma=(1/(n*deltat0))*(f(:,1)'*f(:,1)+f(:,2)'*f(:,2));
end
