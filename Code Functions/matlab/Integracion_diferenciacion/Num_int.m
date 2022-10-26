format compact

fun = @(x) exp(-x.^2).*log(x).^2;

sol = integral(fun,0,Inf)
sol2 = integral(fun,1,5)