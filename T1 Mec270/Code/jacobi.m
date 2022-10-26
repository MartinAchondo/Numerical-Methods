

function [xf,it,rx] = jacobi(A,b,x0,es)

    n = length(b);

    iter = 0;
    x = x0;
    imax = 1000;

    while iter<imax
        x_old = x;
        s = 0;
        for j= 2:n
            s = s + A(1,j)*x(j);
        end
        x(1) = (b(1) - s)/A(1,1);

        for i=2:n
            s = 0;
            for k=1:n
                if k~=i
                    s = s + A(i,k)*x(k);
                end
            end
            x(i) = (b(i)-s)/A(i,i);
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