 [x,y] = meshgrid(-2:.2:2,-2:.2:2) ; 
 
 VectorX = [x./(2*pi*(x.^2+y.^2))]; 
 VectorY = [y./(2*pi*(x.^2+y.^2))]; 
 
 quiver(x,y,VectorX,VectorY)
 axis tight