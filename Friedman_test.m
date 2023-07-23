function s = Friedman_test (n, k, Acc)
a = (12*n/(k*(k+1)))*(sum(Acc.^2)-k*(k+1)^2/4);
b = ((n-1)*a)/(n*(k-1)-a) ;
c = 2.724*sqrt(k*(k+1)/(6*n));
s = [a,b,c];
end
