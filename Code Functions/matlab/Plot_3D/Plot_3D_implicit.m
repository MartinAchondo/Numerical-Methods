format compact

f = @(x,y,z) x.^2 + y.^2 - z.^2;
interval = [-5 5 -5 5 -0 5];

fimplicit3(f,interval)