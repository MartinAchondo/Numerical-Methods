

function [x,iter] = broyden_method(F,J,x0,linear_method,tol)

    imax = 1000;
    xk_old = x0;
    ea = 1;
    iter = 0;
    Bk = J(x0);
    while ea > tol && iter < imax
        sk = linear_method(Bk,-F(xk_old));
        xk = xk_old + sk;
        yk = F(xk) - F(xk_old);
        M = (yk - Bk*sk)*sk';
        Bk = Bk + M/dot(sk,sk);

        ea = norm(sk)/norm(xk);

        xk_old = xk;
        iter = iter + 1;
        
    end
    x = real(xk);
    %x = xk;
end