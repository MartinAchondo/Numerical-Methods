
format compact

tic
lf1 = all_lagrange(0,3*pi,@f1);
toc

tic
nf1 = all_newton(0,3*pi,@f1);
toc

tic
lf2 = all_lagrange(-5,5,@f2);
toc

tic
nf2 = all_newton(-5,5,@f2);
toc

tic
all_lineal(0,3*pi,@f1,[4,10,15,34]);
toc

tic
all_lineal(-5,5,@f2,[4,10,15,35]);
toc

function [y] = f1(x)
    y = sin(x);
end

function [y] = f2(x)
    y = 1./(1+x.^2);
end


function [Y_Lagrange] = all_lagrange(a,b,f)

    pts = 100;
    Y_Lagrange = zeros(7,pts);
    for i=1:7
        order = i + 1;
        n = order + 1;
        x = linspace(a,b,n);
        y = f(x);
        [y_inter,~] = inter_lagrange(x,y,pts);
        Y_Lagrange(n-2,:) = y_inter;
    end
end

function [Y_Newton] = all_newton(a,b,f)

    pts = 100;
    Y_Newton = zeros(7,pts);
    for i=1:7
        order = i + 1;
        n = order + 1;
        x = linspace(a,b,n);
        y = f(x);
        [y_inter,~] = inter_newton(x,y,pts);
        Y_Newton(n-2,:) = y_inter;
    end
end

function [Y_Lineal] = all_lineal(a,b,f,nodes)

    pts = 100;
    Y_Lineal = zeros(7,pts);
    for r=1:4
        i = nodes(r);
        n = i;
        x = linspace(a,b,n);
        y = f(x);
        [y_inter,~] = inter_lineal(x,y,pts);
        Y_Lineal(r,:) = y_inter; 
    end
end