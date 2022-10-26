
format compact


plot_lineal(0,3*pi,@f1,1,'f(x)=sin(x)',[6 12 20 34])

plot_lineal(-5,5,@f2,5,'f(x)=1/(1+x^2)',[5 11 21 35])

plot_lineal_error(0,3*pi,@f1,9,'sin(x)',[6 12 20 34])

plot_lineal_error(-5,5,@f2,13,'f(x)=1/(1+x^2)',[5 11 21 35])


function [y] = f1(x)
    y = sin(x);
end

function [y] = f2(x)
    y = 1./(1+x.^2);
end


function [flag] = plot_lineal(a,b,f,pp,lab,nodes)

    pts = 500;
    x_lineal = linspace(a,b,pts);
    Y_Lineal = zeros(7,pts);
    j = pp;
    for r=1:4
        i = nodes(r);
        n = i;
        x = linspace(a,b,n);
        y = f(x);
        [y_inter,~] = inter_lineal(x,y,pts);
        Y_Lineal(r,:) = y_inter;
    
        x_real = x_lineal;
        figure(j)
        plot(x_real,f(x_real),'k--','DisplayName',lab)
        hold on
        x2 = linspace(a,b,n);
        plot(x_real,Y_Lineal(r,:),'r','DisplayName',strcat('Pol Nodes ',int2str(n)))
        hold on
        xlim([a b])
        xlabel('x')
        ylabel('y')
        scatter(x2,f(x2),'r','DisplayName',strcat('Nodes  ',int2str(n)))
        hold on
        legend
        j = j + 1;    
    end
    flag = true;

end

function [flag] = plot_lineal_error(a,b,f,pp,lab,nodes)

    pts = 500;
    x_lineal = linspace(a,b,pts);
    Y_Lineal = zeros(7,pts);
    j = pp;
    for r=1:4
        i = nodes(r);
        n = i;
        x = linspace(a,b,n);
        y = f(x);
        [y_inter,~] = inter_lineal(x,y,pts);
        Y_Lineal(r,:) = y_inter;

        x_real = x_lineal;
        figure(j)
        x2 = linspace(a,b,n);
        err = f(x_real)-Y_Lineal(r,:);
        disp(strcat('Error of Pol Nodes',int2str(n)))
        e = max(abs(err))
        plot(x_real,err,'r','DisplayName',strcat('Pol Nodes ',int2str(n)))
        hold on
        xlim([a b])
        xlabel('x')
        ylabel('y')
        scatter(x2,f(x2)*0,'r','DisplayName',strcat('Nodes  ',int2str(n)))
        hold on
        legend
        j = j + 1;
      

    end
    flag = true;

end
