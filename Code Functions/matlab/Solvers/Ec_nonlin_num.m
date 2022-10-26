format compact 
syms a 

eq1 = 2*a^2-exp(4)==3;

x = [a];
%se especifica las variables al final
[sol_a] = solve([eq1],x);

a_b = double([sol_a]);
display(a_b)


