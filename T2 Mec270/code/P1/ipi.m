%IPI

function [lambda,w] = ipi(A,linear_method)
    maxiter=100;
    tol=0.0001;
    m=0;
    mu=1;
    x_inicial=[1;2;3];
    w=x_inicial/norm(x_inicial);
    lambda=(transpose(w))*A*w;
    while (norm(A*w-lambda*w))>tol && (m<maxiter)
        Iden=diag([1 1 1]);
        v=linear_method(A+mu*Iden,w);
        w=v/norm(v);
        lambda=(transpose(w))*A*w;
        m=m+1;
        if (norm(A*w-lambda*w))>tol
            m;
        end
    end
end


