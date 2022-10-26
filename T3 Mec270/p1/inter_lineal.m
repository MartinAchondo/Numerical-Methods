
function [y_inter,x_inter] = inter_lineal(x,y,pts)

    n = length(x);
    x_inter = linspace(min(x),max(x),pts);
    y_inter = zeros(1,pts);

    [m,b] = splines_lin(x,y,n);

    i = 1;
    for k=1:pts
        if(x_inter(k)<=x(i+1))
            y_inter(k) = b(i) + m(i)*(x_inter(k)-x(i));
        else
            y_inter(k) = b(i+1) + m(i+1)*(x_inter(k)-x(i+1));
            i = i + 1;
        end
    end

    function [mr,br] = splines_lin(x,y,n)

        mr = zeros(1,n);
        br = zeros(1,n);

        for l=1:n-1
            mr(l) = (y(l+1)-y(l))/(x(l+1)-x(l));
            br(l) = y(l);
        end
    end
end