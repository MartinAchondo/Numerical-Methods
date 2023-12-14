

function [x,iter] = conjgrad_method(F,J,Cf,x0,tol)

    imax = 1000;
    xk_old = x0;
    ea = 1;
    iter = 0;
    rk_old = J(xk_old)'*F(x0);
    pk_old = -rk_old;
    C = Cf();
    while ea > tol && iter < imax

        A = zeros([10 10]);
        for g=1:10
            A(g,g) = dot(C(:,g),F(xk_old));
        end
        A = A + J(xk_old)'*J(xk_old);

        ak = -dot(rk_old,pk_old)/(pk_old'*A*pk_old);
        xk = xk_old + ak*pk_old;
        rk = J(xk)'*F(xk);
        yk = rk - rk_old;
        bk = dot(rk,yk)/dot(rk_old,rk_old);
        pk = -rk + bk*pk_old;

        pk_old = pk;
        rk_old = rk;
        ea = norm(xk-xk_old)/norm(xk);

        xk_old = xk;
        iter = iter + 1;

        if mod(iter,10) == 0
            pk_old = J(xk_old)'*F(xk_old);
        end
    end
    x = xk;
end