
function [y_inter,x_inter] = inter_lagrange(x,y,pts)

    n = length(x);
    x_inter = linspace(min(x),max(x),pts);
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