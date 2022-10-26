 [x,y,z] = meshgrid(-2:.2:2,-2:.2:2,-2:.2:2) ; 
 
 VectorX = x; 
 VectorY = 2.*y; 
 VectorZ = -z; 
 
 quiver3(x,y,z,VectorX,VectorY,VectorZ)
 axis tight