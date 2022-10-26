format compact

x = linspace(0,100,100);

for n=1:100
    y = n.*x;
    
    plot(x,y);
    axis([0 100 0 100]);
    pause(0.1);
end
