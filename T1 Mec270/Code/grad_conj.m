

function [xf,it,rx] = grad_conj(A,b,x0,es)

    iter = 0;

    x = x0;
    s = b - A*x0;
    r = transpose(A)*s;
    p = r;
    q = A*p;

    imax = 1000;
   

    while iter<imax
        x_old = x;
        r_old = r;
        alpha = dot(r_old,r_old)/dot(q,q);
        x = x_old + alpha*p;
        s = s - alpha*q;
        r = transpose(A)*s;
        beta = dot(r,r)/dot(r_old,r_old);
        p = r + beta*p;
        q = A*p;

        normr = norm(x-x_old);

        if normr<es 
            break
        end

        iter = iter + 1;
    end 
    xf = x;
    it = iter;
    rx = b - A*x;
end