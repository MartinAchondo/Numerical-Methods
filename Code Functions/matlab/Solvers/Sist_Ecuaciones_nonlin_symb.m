format compact 
syms a b e

eq1 = 2*a+exp(b)==3+e;
eq2 = 23*a^2 ==2;
x = [a,b];
%se especifica las variables al final
[sol_a,sol_b] = solve([eq1,eq2],x);


display(sol_a);
display(sol_b);


