
function [y_inter,x_inter,y_chebs,x_chebs] = inter_chebs(x,pts,f_eval)

    n = length(x) - 1;
    x_chebs = nodos_chebs(min(x),max(x),n);
    y_chebs = f_eval(x_chebs);
    [y_inter,x_inter] = interpol_lagrange(x_chebs,y_chebs,pts,min(x),max(x));

    function [x_n] = nodos_chebs(a,b,n)
        x_n = zeros(1,n+1);
        for i=0:n
            x_n(i+1) = (a+b)/2 + (b-a)/2*cos(pi*(2*(n-i)+1)/(2*n+2));
        end
    end 

    function [y_inter,x_inter] = interpol_lagrange(x,y,pts,a,b)

        n = length(x);
        x_inter = linspace(a,b,pts);
        y_inter = zeros(1,pts);
        for k=1:pts
            P = 0;
            for i=1:n
                L = L_lagrange(x_inter(k),x,n,i);
                P = P + L*y(i);
            end
            y_inter(k) = P;
        end
    
        function [L] = L_lagrange(x_eval,x_data,n,i)
            L = 1;
            if(x_eval==x_data(i))
                L = 1;
            else
                for j=1:n
                    if (j ~= i)
                        L = L * (x_eval-x_data(j))/(x_data(i)-x_data(j));
                    end
                end
            end
        end
    end
end