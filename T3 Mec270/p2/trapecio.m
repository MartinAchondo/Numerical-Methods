%Metodo del trapecio
function [valor3] = trapecio(n,x,f,k,h)
    suma=0;
    for i=1:n
        y1=eval(subs(f,k,x(i)));
        y2=eval(subs(f,k,x(i+1)));
        suma=suma + h*((y2+y1)/2);
    end
    valor3=suma;
end