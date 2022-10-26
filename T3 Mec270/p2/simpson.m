%Metodo de simpson 
function [valor4] = simpson(n,x,f,k,h)
    v=n;
    suma=0;
    dat1=n/2;
    dat2=floor(n/2);
    if (dat1-dat2)>0 & n>1;
        f1=eval(subs(f,k,x(n-2)));
        f2=eval(subs(f,k,x(n-1)));
        f3=eval(subs(f,k,x(n)));
        f4=eval(subs(f,k,x(n+1)));
        [val1]= simpson38(h,f1,f2,f3,f4);
        suma = suma + val1;
        v = n-3;
    end
    if v >1;
        [val2]= simpson13(h,v,f,k,x);
        suma = suma + val2;
    end
    valor4=suma;

    
end