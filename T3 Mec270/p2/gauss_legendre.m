
%Metodo de gauss-legendre
function [valor6] = gauss_legendre(n,f,k,h,a,b)
    suma=0;
    valor6 = 0;
    if n<7
        [z,w]=gauleg(n);
        for i=1:length(z)
            y=(b+a)/2 + (b-a)*z(i)/2;
            suma=suma+ w(i)*(eval(subs(f,k,y)));
        end
        valor6 =((b-a)/2)*suma;
    end
    
    %Calculo de w y z para gauss-legendre
    function [z,w] = gauleg(n)
    if n==1
        x= [0] ;
        c=[2];
    end
    if n==2
        x=[-0.577350269 0.577350269];
        c=[1 1];
    end
    if n==3
        x=[-0.774596669 0.0 0.774596669]; 
        c=[0.5555556 0.8888889 0.5555556];
    end
    if n==4
        x=[-0.861136312 -0.339981044 0.339981044 0.861136312];
        c=[0.3478548 0.6521452 0.6521452 0.3478548];
    end
    if n==5
        x=[-0.906179846 -0.538469310 0.0 0.538469310 0.906179846]; 
        c=[0.2369269 0.4786287 0.5688889 0.4786287 0.2369269];
    end
    if n==6
        x=[-0.932469514 -0.661209386 -0.238619186 0.238619186 0.661209386 0.932469514];
        c=[0.1713245 0.3607616 0.4679139 0.4679139 0.3607616 0.1713245];
    end
    z=x;
    w=c;
end

end