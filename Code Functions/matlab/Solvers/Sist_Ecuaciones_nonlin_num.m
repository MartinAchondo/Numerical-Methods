format compact 
syms a b 

eq1 = 2*a+exp(b)==3;
eq2 = 23*a^2 ==2;
x = [a,b];
%se especifica las variables al final
[sol_a,sol_b] = solve([eq1,eq2],x);

a_b = double([sol_a,sol_b]);

for i=1:length(a_b)
    display("Sol" + i + " = " )
    display(a_b(i))
    display(a_b(i+length(x)))
end


