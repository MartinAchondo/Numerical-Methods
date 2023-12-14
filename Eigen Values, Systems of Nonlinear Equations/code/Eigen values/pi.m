
function [lambda,w] = pi(A)
    maxiter=100;
    tol=0.0001;
    m=0;
    x_inicial=[1;2;3];
    w=x_inicial/norm(x_inicial);
    lambda=(transpose(w))*A*w;
    while (norm(A*w-lambda*w))>tol && (m<maxiter)
        w=A*w;
        w=w/norm(w);
        lambda=(transpose(w))*A*w;
        m=m+1;
    end
end



