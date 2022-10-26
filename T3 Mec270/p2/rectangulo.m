
%Metodo del rectangulo izquierda
function [valor1] = rectangulo(n,x,f,k,h)
    suma=0;
    for i=1:n
        y=eval(subs(f,k,x(i)));
        suma=suma + h*y;
    end
    valor1=suma;
end