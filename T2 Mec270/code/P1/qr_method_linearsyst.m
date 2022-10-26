
function [x] = qr_method_linearsyst(A,b)

    [Qt,R] = qr_factor(A);
    b = Qt*b;
    n = size(b);
    x = zeros([n 1]);
    for i=n:-1:1
        x(i) = b(i);
        for j=i+1:n
            x(i) = x(i) - R(i,j)*x(j);
        end
        x(i) = 1/R(i,i) * x(i);
    end

    function [Qt,R] = qr_factor(A)
        [m,n] = size(A);
        H = eye(n);
        Hx = H;
        for k=1:n-1
            u = A(k:n,k);
            u(1) = u(1) + sign(A(k,k))*norm(A(k:n,k));
            v = u/u(1);
            beta = 2/(v'*v);
            Hk = eye(n+1-k) - beta*v*v';
            H = eye(n);
            H(k:m,k:m) = Hk;
            A = H*A;
            Hx = H*Hx;
        end
        R = A;
        Qt = Hx;
    end
end



