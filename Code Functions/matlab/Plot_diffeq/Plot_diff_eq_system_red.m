[x,y]=meshgrid(-5:.4:5,-10:.51:10);

dy=-9.81*sin(x).*sign(y).^2;
dx = y.*sign(x).^2;

dxu=dx./sqrt(dx.^2+dy.^2);
dyu=dy./sqrt(dx.^2+dy.^2);

quiver(x,y,dxu,dyu)
axis tight