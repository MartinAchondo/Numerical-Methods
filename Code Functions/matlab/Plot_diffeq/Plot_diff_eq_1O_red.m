[x,y]=meshgrid(-3:.3:3,-2:.3:2);
dy = x+sin(y);
dx=ones(size(dy));

quiver(x,y,dx,dy)
axis tight