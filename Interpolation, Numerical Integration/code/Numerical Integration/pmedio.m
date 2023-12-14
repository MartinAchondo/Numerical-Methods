
%Metodo del punto medio
function [valor2] = pmedio(n,x,f,k,h)
    suma=0;
    for i=1:n;
        j=(x(i)+x(i+1))/2;
        y=eval(subs(f,k,j));
        suma=suma + h*y;
    end
    valor2=suma;
end