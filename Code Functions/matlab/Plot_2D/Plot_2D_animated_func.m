format compact

T1 = 45;
T0 = 12;
theta_0 = T0-T1;
alpha = 7.3*10^(-5);
l = 1;

x = linspace(0,1,100);
n=200;
for t=1:200
    y = fr(n,T1,theta_0,alpha,t,l,x);
    
    plot(x,y);
    pause(0.1);
end

function res = fr(n,T1,theta_0,alpha,t,l,x)
    theta = zeros(1,100);
    theta = theta + T1;
    for k = 1:n
        theta = theta + 4/((2*k-1)*pi)*theta_0*exp(-alpha*(k*pi/l)^2*t)*sin((2*k-1)*pi.*x/l); 
    end
    res = theta; 
end
