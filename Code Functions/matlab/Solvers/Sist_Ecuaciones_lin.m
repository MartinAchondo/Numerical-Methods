syms x1 x2 x3 a b ;
%sistema lineal!!
eq1 = x1 + 2*x2 +a*b^3 == a;
eq2 = x3+25*x2 == b^2;
eq3 = x1-x2 + x3 == a;

eq = [eq1,eq2,eq3];
x = [x1,x2,x3];

[A,B] = equationsToMatrix(eq,x);
sol = linsolve(A,B);

for i=1:length(x)
   display("x"+i+" = "+ string(sol(i)));  
end