

function [vp,VP] = qr_hessenberg(A)

    iter=1000;
    tol = 10^-6;
    [m,n]=size(A);
    Vp=eye(n);

    %calculo matriz de hessenberg
    Q=eye(m);
    for i=1:m-2
        H=eye(m);
        x=A(i+1:m,i);
        w=[-sign(x(1))*norm(x);zeros(m-i-1,1)];
        v=w-x;
        l=length(x);
        s=x(2:l)'*x(2:l);
        if (s==0)
            variable=0;
        else
            variable=2/(norm(v)^2);
        end
        H(i+1:m,i+1:m)=eye(i+1)-variable*v*v';
        A=H*A*H;
    end

    %iteracion qr
    T=A;
    T_old = T;
    i = 0;
    error = 1;
    while i<iter && error > tol
        [Q,R]=qr_factor(T_old);
        T=R*Q;
        Vp=Vp*Q;
        error = norm(diag(T_old)-diag(T));
        T_old = T;
        i = i+1;
    end

    vp = T; 
    VP = Vp*Q;
    
    function [Qt,R] = qr_factor(T)
        [m,n] = size(T);
        H = eye(n);
        Hx = H;
        for k=1:n-1
            u = T(k:n,k);
            u(1) = u(1) + sign(T(k,k))*norm(T(k:n,k));
            v = u/u(1);
            beta = 2/(v'*v);
            Hk = eye(n+1-k) - beta*v*v';
            H = eye(n);
            H(k:m,k:m) = Hk;
            T = H*T;
            Hx = H*Hx;
        end
        R = T;
        Qt = Hx;
        Qt=Qt';
    end

end