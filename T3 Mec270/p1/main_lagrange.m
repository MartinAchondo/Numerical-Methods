
format compact


plot_lagrange(0,3*pi,@f1,1,'f(x)=sin(x)')

plot_lagrange(-5,5,@f2,5,'f(x)=1/(1+x^2)')

plot_lagrange_error(0,3*pi,@f1,9,'sin(x)')

plot_lagrange_error(-5,5,@f2,13,'f(x)=1/(1+x^2)')

function [y] = f1(x)
    y = sin(x);
end

function [y] = f2(x)
    y = 1./(1+x.^2);
end


function [flag] = plot_lagrange(a,b,f,pp,lab)

    pts = 100;
    x_Lagrange = linspace(a,b,pts);
    Y_Lagrange = zeros(7,pts);
    j = pp;
    for i=1:7
        order = i + 1;
        n = order + 1;
        x = linspace(a,b,n);
        y = f(x);
        [y_inter,~] = inter_lagrange(x,y,pts);
        Y_Lagrange(n-2,:) = y_inter;
    
        if (mod(order,2)~=0)
            x_real = x_Lagrange;
            figure(j)
            plot(x_real,f(x_real),'k--','DisplayName',lab)
            hold on
            plot(x_real,Y_Lagrange(i-1,:),'b','DisplayName',strcat('Pol Order ',int2str(order-1)))
            hold on
            x1 = linspace(a,b,n-1);
            x2 = linspace(a,b,n);
            scatter(x1,f(x1),'b','DisplayName',strcat('Nodes Order ',int2str(order-1)))
            hold on
            plot(x_real,Y_Lagrange(i,:),'r','DisplayName',strcat('Pol Order ',int2str(order)))
            hold on
            xlim([a b])
            xlabel('x')
            ylabel('y')
            scatter(x2,f(x2),'r','DisplayName',strcat('Nodes Order ',int2str(order)))
            hold on
            legend
            j = j + 1;
        end
    
        if(i==7)
            figure(j)
            plot(x_real,f(x_real),'k--','DisplayName',lab)
            hold on
            plot(x_real,Y_Lagrange(i,:),'b','DisplayName',strcat('Pol Order ',int2str(order)))
            hold on
            xlim([a b])
            xlabel('x')
            ylabel('y')
            x1 = linspace(a,b,n);
            scatter(x1,f(x1),'b','DisplayName',strcat('Nodes Order ',int2str(order)))
            hold on
            legend
        end
    
    end
    flag = true;

end

function [flag] = plot_lagrange_error(a,b,f,pp,lab)

    pts = 100;
    x_Lagrange = linspace(a,b,pts);
    Y_Lagrange = zeros(7,pts);
    j = pp;
    for i=1:7
        order = i + 1;
        n = order + 1;
        x = linspace(a,b,n);
        y = f(x);
        [y_inter,~] = inter_lagrange(x,y,pts);
        Y_Lagrange(n-2,:) = y_inter;
    
        if (mod(order,2)~=0)
            x_real = x_Lagrange;
            figure(j)
            err = f(x_real)-Y_Lagrange(i-1,:);
            disp(strcat('Error of Pol Order',int2str(order-1)))
            e = max(abs(err))
            plot(x_real,err,'b','DisplayName',strcat('Err Pol Order ',int2str(order-1)))
            hold on
            x1 = linspace(a,b,n-1);
            x2 = linspace(a,b,n);
            scatter(x1,f(x1)*0,'b','DisplayName',strcat('Nodes Order ',int2str(order-1)))
            hold on
            err = f(x_real)-Y_Lagrange(i,:);
            disp(strcat('Error of Pol Order',int2str(order)))
            e = max(abs(err))
            plot(x_real,err,'r','DisplayName',strcat('Err Pol Order ',int2str(order)))
            hold on
            xlim([a b])
            xlabel('x')
            ylabel('y')
            scatter(x2,f(x2)*0,'r','DisplayName',strcat('Nodes Order ',int2str(order)))
            hold on
            legend
            j = j + 1;
        end
    
        if(i==7)
            figure(j)
            err = f(x_real)-Y_Lagrange(i,:);
            disp(strcat('Error of Pol Order',int2str(order)))
            e = max(abs(err))
            plot(x_real,err,'b','DisplayName',strcat('Err Pol Order ',int2str(order)))
            hold on
            xlim([a b])
            xlabel('x')
            ylabel('y')
            x1 = linspace(a,b,n);
            scatter(x1,f(x1)*0,'b','DisplayName',strcat('Nodes Order ',int2str(order)))
            hold on
            legend
        end
    
    end
    flag = true;

end
