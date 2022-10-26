%solo cambiar extremos
[X,Y] = meshgrid(-4:.5:4);
Z = Y.*sin(X) - X.*cos(Y);

surf(X,Y,Z,'FaceAlpha',0.5);