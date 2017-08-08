function nb = fct_nb_peak(s)
% Compute the number of extreme peaks in the signal s
%

sstd = std(s,[],1);
iii = bsxfun(@ge, s, 10*sstd);
jjj = (diff(iii,1,1)==1);
nb = sum(jjj,1);