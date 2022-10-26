clc
clear all
tic
a=0;
b=1;
syms k;
f= [exp(-k^{2})];
g= [-2*k*(exp(-k^{2}))];
I=int(f,a,b);
n=6;
h=(b-a)/n;
x=a:h:b;

%calculo de los metodos 
[rectang] = rectangulo(n,x,f,k,h) 
[pmedioo] = pmedio(n,x,f,k,h)
[trap] = trapecio(n,x,f,k,h)
[simps] = simpson(n,x,f,k,h)
[trap_co]=trapecio_corregido(n,x,f,g,k,h,a,b)
[Gauss_leg]= gauss_legendre(n,f,k,h,a,b)

%Calculo de errores
RI=eval(I);
error_Rectangulo= abs(RI-rectang)
error_Pmedio= abs(RI-pmedioo)
error_trapecio= abs(RI-trap)
error_simpson= abs(RI-simps)
error_trapecio_Corregido = abs(RI-trap_co)
error_Gauss_legendre= abs(RI-Gauss_leg)
toc
