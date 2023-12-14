

function [x] = cholesky_method(A,b)

    [L] = cholesky_fact(A);

    y = fow_subs(L,b);
    x = back_subs(L',y);

    function [x] = fow_subs(L,b)
        [~,n] = size(L);
        x = zeros(n,1);
        for i=1:n
            x(i) = b(i);
            for j=1:i-1
                x(i) = x(i) - L(i,j)*x(j);
            end
            x(i) = 1/L(i,i) * x(i);
        end
    end

    function [x] = back_subs(U,b)
        [~,n] = size(U);
        x = zeros(n,1);
        for i=n:-1:1
            x(i) = b(i);
            for j=i+1:n
                x(i) = x(i) - U(i,j)*x(j);
            end
            x(i) = 1/U(i,i) * x(i);
        end
    end

    function [L] = cholesky_fact(A)

        [~,n] = size(A);
        L = zeros(3);
    
        L(1,1) = sqrt(A(1,1));
    
        for k=1:n
            for i=1:k-1
                s = 0;
                for j=1:i-1
                    s = s + L(i,j)*L(k,j);
                end
                L(k,i) = (A(k,i) - s)/L(i,i);
            end
            s = 0;
            L(k,k) = sqrt(A(k,k) - norm(L(k,1:k-1))^2);

        end


     
    end

end




