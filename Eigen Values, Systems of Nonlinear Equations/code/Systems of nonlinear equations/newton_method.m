

function [x,iter] = newton_method(F,J,x0,linear_method,tol)

    imax = 1000;
    xk_old = x0;
    ea = 1;
    iter = 0;
    while ea > tol && iter < imax
        sk = linear_method(J(xk_old),-F(xk_old));
        xk = xk_old + sk;
        ea = norm(sk)/norm(xk);

        xk_old = xk;
        iter = iter + 1;
    end
    x = real(xk);
    %x = xk;
end
