
function [y_inter,x_inter] = inter_newton(x,y,pts)
    
    n = length(x);
    lam = fin_dif(x,y);
    x_inter = linspace(min(x),max(x),pts);
    y_inter = zeros(1,pts);

    for z=1:pts
        P = 0;
        for l=1:n
            Nj = N_newton(x_inter(z),x,l);
            P = P + lam(l)*Nj;
        end
        y_inter(z) = P;
    end

    function [N] = N_newton(x_eval,x_data,l)
        N = 1;
        if (l==1)
            N = 1;
        else
            for q=1:l-1
                N = N * (x_eval - x_data(q));
            end
        end

    end


    function [lambda] = fin_dif(x,y)
       
        n = length(x);
        fn = zeros(n,n);
        for j=1:n
            fn(j,j) = y(j); 
        end
    
        for i=2:n
            for k=1:n-i+1
                fn(k,i+k-1) = (fn(k+1,i+k-1)-fn(k,i+k-2))/(x(i+k-1)-x(k));
            end
        end

        lambda = fn(1,:);

    end

end