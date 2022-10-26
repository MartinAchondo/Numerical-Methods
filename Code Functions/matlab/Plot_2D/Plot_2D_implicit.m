f = @(x,y) x.^2 + y.^2 - 25;
interval = [-5 5 -5 5];

fimplicit(f,interval)