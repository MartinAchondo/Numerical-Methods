format compact 
syms a b e

eq1 = 2*a+exp(b)==3+e;

x = [a];
%se especifica las variables al final
[sol_a] = solve([eq1],x);


display(sol_a);