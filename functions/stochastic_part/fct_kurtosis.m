function kurt = fct_kurtosis(s)
% Compute the kurtosis of the signal s
% along the line
%

s = bsxfun(@plus, s, -mean(s,1));
s = bsxfun(@times,s, 1./std(s,[],1));
% s = (s-mean(s,1))/std(s);
kurt = mean(s.^4,1);