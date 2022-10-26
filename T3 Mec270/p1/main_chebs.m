
format compact


plot_tchebyshev(0,3*pi,@f1,1,'f(x)=sin(x)')

plot_tchebyshev(-5,5,@f2,5,'f(x)=1/(1+x^2)')

plot_tchebyshev_error(0,3*pi,@f1,9,'sin(x)')

plot_tchebyshev_error(-5,5,@f2,13,'f(x)=1/(1+x^2)')

function [y] = f1(x)
    y = sin(x);
end

function [y] = f2(x)
    y = 1./(1+x.^2);
end


function [flag] = plot_tchebyshev(a,b,f,pp,lab)

    pts = 100;
    x_chebs = linspace(a,b,pts);
    Y_Chebs = zeros(7,pts);
    j = pp;
    for i=1:7
        order = i + 1;
        n = order + 1;
        x = linspace(a,b,n);
        y = f(x);
        [y_inter,x2,yi2,xi2] = inter_chebs(x,pts,f);
        Y_Chebs(n-2,:) = y_inter;

        if(mod(order,2)==0)
            xi1 = xi2;
            yi1 = yi2;
            x1 = x2;
        end

        if (mod(order,2)~=0)
            x_real = x_chebs;
            figure(j)
            plot(x_real,f(x_real),'k--','DisplayName',lab)
            hold on
            plot(x1,Y_Chebs(i-1,:),'b','DisplayName',strcat('Pol Order ',int2str(order-1)))
            hold on
            scatter(xi1,f(xi1),'b','DisplayName',strcat('Nodes Order ',int2str(order-1)))
            hold on
            plot(x2,Y_Chebs(i,:),'r','DisplayName',strcat('Pol Order ',int2str(order)))
            hold on
            xlim([a b])
            xlabel('x')
            ylabel('y')
            scatter(xi2,f(xi2),'r','DisplayName',strcat('Nodes Order ',int2str(order)))
            hold on
            legend
            j = j + 1;
        end
    
        if(i==7)
            figure(j)
            plot(x1,f(x1),'k--','DisplayName',lab)
            hold on
            plot(x1,Y_Chebs(i,:),'b','DisplayName',strcat('Pol Order ',int2str(order)))
            hold on
            xlim([a b])
            xlabel('x')
            ylabel('y')
            x1 = linspace(a,b,n);
            scatter(xi1,f(xi1),'b','DisplayName',strcat('Nodes Order ',int2str(order)))
            hold on
            legend
        end
    
    end
    flag = true;

end

function [flag] = plot_tchebyshev_error(a,b,f,pp,lab)

    pts = 100;
    x_chebs = linspace(a,b,pts);
    Y_Chebs = zeros(7,pts);
    j = pp;
    for i=1:7
        order = i + 1;
        n = order + 1;
        x = linspace(a,b,n);
        y = f(x);
        [y_inter,x2,yi2,xi2] = inter_chebs(x,pts,f);
        Y_Chebs(n-2,:) = y_inter;

        if(mod(order,2)==0)
            xi1 = xi2;
            yi1 = yi2;
            x1 = x2;
        end
    
        if (mod(order,2)~=0)
            x_real = x_chebs;
            figure(j)
            err = f(x1)-Y_Chebs(i-1,:);
            disp(strcat('Error of Pol Order',int2str(order-1)))
            e = max(abs(err))
            plot(x1,err,'b','DisplayName',strcat('Err Pol Order ',int2str(order-1)))
            hold on
            scatter(xi1,f(xi1)*0,'b','DisplayName',strcat('Nodes Order ',int2str(order-1)))
            hold on
            err = f(x2)-Y_Chebs(i,:);
            disp(strcat('Error of Pol Order',int2str(order)))
            e = max(abs(err))
            plot(x2,err,'r','DisplayName',strcat('Err Pol Order ',int2str(order)))
            hold on
            xlim([a b])
            xlabel('x')
            ylabel('y')
            scatter(xi2,f(xi2)*0,'r','DisplayName',strcat('Nodes Order ',int2str(order)))
            hold on
            legend
            j = j + 1;
        end
    
        if(i==7)
            figure(j)
            err = f(x1)-Y_Chebs(i,:);
            disp(strcat('Error of Pol Order',int2str(order)))
            e = max(abs(err))
            plot(x1,err,'b','DisplayName',strcat('Err Pol Order ',int2str(order)))
            hold on
            xlim([a b])
            xlabel('x')
            ylabel('y')
            scatter(xi1,f(xi1)*0,'b','DisplayName',strcat('Nodes Order ',int2str(order)))
            hold on
            legend
        end
    
    end
    flag = true;

end
