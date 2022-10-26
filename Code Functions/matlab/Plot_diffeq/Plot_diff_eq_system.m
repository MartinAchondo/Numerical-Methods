%dx=dy/dt  dy=dy/dt
[x,y] = meshgrid(-2:0.2:2);
dx = x.^2-3*x.*y+y;
dy = -5*x+sin(x.*y);

r = ( dx.^2 + dy.^2 ).^0.5;
px = dx./r;
py = dy./r;

quiver(x,y,px,py);

axis tight