function corr=correlacion(f,g,mediaf,mediag)

cov=covarianza(f,g,mediaf,mediag);
varf=varianza(f);
varg=varianza(g);
corr=cov/(sqrt(varf*varg));
end