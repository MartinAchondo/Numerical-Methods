function [val2]= simpson13(h,v,f,k,x)
    n=v;
    suma= eval(subs(f,k,x(1)));
    for i=1:2:n-2;
        suma = suma+ 4*eval(subs(f,k,x(i+1)))+ 2*eval(subs(f,k,x(i+2)));
    end
    suma = suma + 4*eval(subs(f,k,x(n))) + eval(subs(f,k,x(n+1)));
    val2 = (h*suma)/3;
end