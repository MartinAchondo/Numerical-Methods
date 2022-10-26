M = 100 ; N = 100 ;
r = linspace(0,100,M) ; 
t = linspace(0,2*pi,N) ;

[R,T] = meshgrid(r,t) ;
X = R.*cos(T) ;
Y = R.*sin(T) ;
Z = R.^3/6 ;

surf(X,Y,Z) ;