
function [x] = lu_method(A,B)
    AA=A;
    n= size(A,2);
    tole=1E-8;
    O=zeros(1,n);
    S=zeros(1,n);
    X=[0;0;0];
    error=0;

    % primer paso
    for i=1:n
        O(i)=i;
        S(i)=abs(A(i,1));
        for j=2:n;
            if abs(A(i,j))>S(i);
                S(i)= abs(A(i,j));
            end
        end
    end

    for k=1:n-1;
        p=k;
        bigg=abs(A(O(k),k)/S(O(k)));
        for a=k+1:n
            dato= abs(A(O(a),k)/S(O(a)));
            if dato> bigg;
                bigg= dato;
                p=a;
            end
        end
        dato=O(p);
        O(p)=O(k);
        O(k)=dato;
        if abs(A(O(k),k)/S(O(k)))<tole;
            error=-1;
            disp(A(O(k),k)/S(O(k)))
            break
        end
        for i = k+1:n
            factor=A(O(i),k)/A(O(k),k);
            A(O(i),k)=factor;
            for j= k+1:n
                A(O(i),j)=A(O(i),j)-factor*A(O(k),j);
            end
        end
    end
    if abs(A(O(k),k)/S(O(k)))<tole;
        error=-1;
        disp(A(O(k),k)/S(O(k)))
    end

    %segundo paso
    if error~=-1
        for i=2:n
            valor1=B(O(i));
            for j=1:i-1;
                valor1=valor1-A(O(i),j)*B(O(j));
            end
            B(O(i))=valor1;
        end
        X(n)= B(O(n))/A(O(n),n);
        for i=n-1:-1:1;
            valor2=0;
            for j = i+1:n
                valor2=valor2+A(O(i),j)*X(j);
            end
            X(i)=(B(O(i))-valor2)/A(O(i),i);
        end
    end

    x = X;

end