function var=varianza(f)
mediaf=esperanza(f);
var=covarianza(f,f,mediaf,mediaf);
end