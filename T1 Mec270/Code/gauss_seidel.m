

function [xf,it,rx] = gauss_seidel(A,b,x0,es)

    n = length(b);
    iter = 0;

    x = x0;
    imax = 1000;

    while iter<imax
        x_old = x;
        for l = 1:n
            sum1 = 0;
            sum2 = 0;
            l1 = l-1;
            if l ~= 1
                for j1 = 1:l1
                    sum1 = sum1 + A(l,j1)*x(j1);
                end
            end
            l2 = l+1;
            if l ~= n
                for j2 = l2:n
                    sum2 = sum2 + A(l,j2)*x_old(j2);
                end
            end
            sum = sum1 + sum2;
            x(l) =  (1/A(l,l))*(b(l) - sum);
        end
        normr = norm(x-x_old);

        if normr<es 
            break
        end
        iter = iter + 1;
    end
    xf = x;
    it = iter;
    rx = b-A*x;
end
