format compact

[y,x] = meshgrid(0:.2:4,0:.2:4);

m = sin(y); 
L = sqrt(1+m.^2);

quiver(x,y,1./L,m./L)

axis tight