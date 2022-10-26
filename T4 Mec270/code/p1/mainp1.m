clc
clear all
Pr=10;
beta=8/3;
valor=0:2:30;
valor1=[0.5,10,28];

t0=0;
tf=100;
n = 10000;
h=(tf-t0)/n;
t = 0:h:n*h-h;


%graficos

%for i=1:15
%    Ra=valor(i);
%    [x,y,z]=Kutta(Ra,Pr,beta,n,h);
%    figure(1)
%    nexttile
%    ploteo=plot(t,z)
%    xlabel('t')
%    ylabel('z')
%    title('z vs t para Ra:',Ra)
%end

%for i=1:3
%    Ra=valor1(i);
%    [x,y,z]=Kutta(Ra,Pr,beta,n,h);
%    figure(i)
%    nexttile
%    ploteo=plot3(x,y,z)
%    comet3(x,y,z)
%    xlabel('x')
%    ylabel('y')
%    zlabel('z')
%    title('Espacio (x,y,z) para Ra:',Ra)
%    grid on
%end

%for i=1:3
%    Ra=valor1(i);
%    [x,y,z]=Kutta(Ra,Pr,beta,n,h);
%    figure(i)
%    nexttile
%    subplot(2,3,1)
%    plot(x,y)
%    xlabel('x')
%    ylabel('y')
%    title('x vs y para Ra:',Ra)
%    grid on
%    subplot(2,3,2)
%    plot(x,z)
%    xlabel('x')
%    ylabel('z')
%    title('x vs z para Ra:',Ra)
%    grid on
%    subplot(2,3,3)
%    plot(y,z)
%    xlabel('y')
%    ylabel('z')
%    title('y vs z para Ra:',Ra)
%    grid on   
%end

figure(6)
[x,y,z]=Kutta(24,Pr,beta,n,h);
comet3(x,y,z)
xlabel('x')
ylabel('y')
zlabel('z')
title('Espacio (x,y,z) para Ra:',24)
grid on


function [x,y,z]=Kutta(Ra,Pr,beta,n,h)
x1=1;
y1=1;
z1=1;

x = zeros(n,1);
y = zeros(n,1);
z = zeros(n,1);

x(1)=x1;
y(1)=y1;
z(1)=z1;

for i=1:n-1
    k1=[fun_dx(Pr,x(i),y(i)) fun_dy(Ra,x(i),y(i),z(i)) fun_dz(beta,y(i),x(i),z(i))];
    k2=[fun_dx(Pr,x(i)+(h/2)*k1(1),y(i)+(h/2)*k1(2)) fun_dy(Ra,x(i)+(h/2)*k1(1),y(i)+(h/2)*k1(2),z(i)+(h/2)*k1(3)) fun_dz(beta,y(i)+(h/2)*k1(2),x(i)+(h/2)*k1(1),z(i)+(h/2)*k1(3))];
    k3=[fun_dx(Pr,x(i)+(h/2)*k2(1),y(i)+(h/2)*k2(2)) fun_dy(Ra,x(i)+(h/2)*k2(1),y(i)+(h/2)*k2(2),z(i)+(h/2)*k2(3)) fun_dz(beta,y(i)+(h/2)*k2(2),x(i)+(h/2)*k2(1),z(i)+(h/2)*k2(3))];
    k4=[fun_dx(Pr,x(i)+h*k3(1),y(i)+h*k3(2)) fun_dy(Ra,x(i)+h*k3(1),y(i)+h*k3(2),z(i)+h*k3(3)) fun_dz(beta,y(i)+h*k3(2),x(i)+h*k3(1),z(i)+(h/2)*k3(3))];
    x(i+1)=x(i)+(h/6)*(k1(1)+2*k2(1)+2*k3(1)+k4(1));
    y(i+1)=y(i)+(h/6)*(k1(2)+2*k2(2)+2*k3(2)+k4(2));
    z(i+1)=z(i)+(h/6)*(k1(3)+2*k2(3)+2*k3(3)+k4(3));
end

x;
y;
z;
end

%funciones derivadas
function [u]=fun_dx(Pr,x,y)
    u=Pr*(y-x);
end

function [v]=fun_dy(Ra,x,y,z)
    v=Ra*x-y-x*z;
end
function [w]=fun_dz(beta,y,x,z)
    w=y*x-beta*z;
end






