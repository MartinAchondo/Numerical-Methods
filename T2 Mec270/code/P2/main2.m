format compact
clear
clc

tol = 10^-7;
x0 = [ 1 ;1 ;1 ;1 ;1 ;1 ;1 ;1 ;1 ;1]*0.5;

options = optimoptions('fsolve','Display','none');
[xr1] = fsolve(@F_hw,x0,options);

disp("matlab H-W")
disp(xr1)
disp("-------")

options = optimoptions('fsolve','Display','none');
[xr2] = fsolve(@F_dw,x0,options);

disp("matlab D-W")
disp(xr2)
disp("-------")

disp("newton")

tic
[xr_n,iter_n] = newton_method(@F_hw,@J_hw,x0,@qr_method_linearsyst,tol);
toc

disp(xr_n)
disp(iter_n)
%disp(F_dw(xr_n))
disp("-------")


disp("broyden")
tic
[xr_b,iter_b] = broyden_method(@F_hw,@J_hw,x0,@qr_method_linearsyst,tol);
toc

disp(xr_b)
disp(iter_b)
%disp(F_hw(xr_b))
disp("-------")

disp("conj grad")
tic
[xr_cg,iter_cg] = conjgrad_method(@F_dw,@J_dw,@dJ_dw,x0,tol);
toc
disp(xr_cg)
disp(iter_cg)
%disp(J_dw(xr_cg))
disp("-------")

A = zeros(10,5);
A(:,1) = xr1;
A(:,2) = xr2;
A(:,3) = xr_n;
A(:,4) = xr_b;
A(:,5) = xr_cg;
Results_xr = real(A);
Results_xr


E_R = zeros(3,1);
for i=3:5
    if i<5
        E_R(i-2) = norm(xr1-A(:,i))/norm(xr1);
    else
        E_R(i-2) = norm(xr2-A(:,i))/norm(xr2);
    end
end

Res = zeros(3,1);
for i=3:5
    if i<5
        Res(i-2) = norm(F_hw(A(:,i))); 
    else
        Res(i-2) = norm(F_dw(A(:,i)));
    end
end

disp("Error")
E_R
disp("Residual")
Res


% funciones para F y J respecto a HW y DW
function [F] = F_hw(x)
    a = 1.852;
    k12 = 351.36;
    k23 = 175.68;
    k34 = 2329.37;
    k45 = k34;
    k56 = 1017.31;
    k67 = 508.66;
    k78 = k67;
    k15 = k23;
    k47 = k56;
    k38 = k47;
    F = [x(1)-x(2)-0.038;
         x(2)-x(3)-x(10)-0.008;
         x(3)-x(9)+x(4)-0.02;
         x(8)-x(4)-x(5)-0.02;
         x(5)-x(6)-0.006;
         x(6)+x(9)-x(7)-0.008;
         x(7)+x(10)-0.1;
         -k12*x(1)^a+k15*x(8)^a+k45*x(4)^a-k34*x(3)^a-k23*x(2)^a;
         -k45*x(4)^a+k56*x(5)^a+k67*x(6)^a-k47*x(9)^a;
         k34*x(3)^a+k47*x(9)^a+k78*x(7)^a-k38*x(10)^a
        ];
end

function [J] = J_hw(x)
    a = 1.852;
    b = a-1;
    k12 = 351.36*a;
    k23 = 175.68*a;
    k34 = 2329.37*a;
    k45 = k34;
    k56 = 1017.31*a;
    k67 = 508.66*a;
    k78 = k67;
    k15 = k23;
    k47 = k56;
    k38 = k47;
    J = [ 1 -1 0 0 0 0 0 0 0 0 ; 
        0 1 -1 0 0 0 0 0 0 -1 ;
        0 0 1 1 0 0 0 0 -1 0;
        0 0 0 -1 -1 0 0 1 0 0 ;
        0 0 0 0 1 -1 0 0 0 0;
        0 0 0 0 0 1 -1 0 1 0;
        0 0 0 0 0 0 1 0 0 1;
        -k12*x(1)^b -k23*x(2)^b -k34*x(3)^b k45*x(4)^b 0 0 0 k15*x(8)^b 0 0;
        0 0 0 -k45*x(4)^b k56*x(5)^b k67*x(6)^b 0 0 -k47*x(9)^b 0;
        0 0 k34*x(3)^b 0 0 0 k78*x(7)^b 0 k47*x(9)^b -k38*x(10)^b
        ];
end


function [F] = F_dw(x)
    a = 2;
    g = 9.81;
    w = 8/(g*pi^2);
    k1 = 0.0215*w*300/(0.255^5);
    k2 = 0.0230*w*150/(0.255^5);
    k3 = 0.0287*w*150/(0.15^5);
    k4 = 0.0237*w*150/(0.15^5);
    k5 = 0.0234*w*300/(0.205^5);
    k6 = 0.0239*w*150/(0.205^5);
    k7 = 0.0229*w*150/(0.205^5);
    k8 = 0.0215*w*150/(0.255^5);
    k9 = 0.0263*w*300/(0.205^5);
    k10 = 0.0234*w*300/(0.205^5);

    F = [x(1)-x(2)-0.038;
         x(2)-x(3)-x(10)-0.008;
         x(3)-x(9)+x(4)-0.02;
         x(8)-x(4)-x(5)-0.02;
         x(5)-x(6)-0.006;
         x(6)+x(9)-x(7)-0.008;
         x(7)+x(10)-0.1;
         -k1*x(1)^a+k8*x(8)^a+k4*x(4)^a-k3*x(3)^a-k2*x(2)^a;
         -k4*x(4)^a+k5*x(5)^a+k6*x(6)^a-k9*x(9)^a;
         k3*x(3)^a+k9*x(9)^a+k7*x(7)^a-k10*x(10)^a
        ];
        F(8) = F(8)/k1;
        F(9) = F(9)/k1;
        F(10) = F(10)/k1;
end

function [J] = J_dw(x)
    a = 2;
    g = 9.81;
    w = 8/(g*pi^2);
    k1 = 0.0215*w*300/(0.255^5)*a;
    k2 = 0.0230*w*150/(0.255^5)*a;
    k3 = 0.0287*w*150/(0.15^5)*a;
    k4 = 0.0237*w*150/(0.15^5)*a;
    k5 = 0.0234*w*300/(0.205^5)*a;
    k6 = 0.0239*w*150/(0.205^5)*a;
    k7 = 0.0229*w*150/(0.205^5)*a;
    k8 = 0.0215*w*150/(0.255^5)*a;
    k9 = 0.0263*w*300/(0.205^5)*a;
    k10 = 0.0234*w*300/(0.205^5)*a;
    J = [   1 -1 0 0 0 0 0 0 0 0 ; 
            0 1 -1 0 0 0 0 0 0 -1 ;
            0 0 1 1 0 0 0 0 -1 0;
            0 0 0 -1 -1 0 0 1 0 0 ;
            0 0 0 0 1 -1 0 0 0 0;
            0 0 0 0 0 1 -1 0 1 0;
            0 0 0 0 0 0 1 0 0 1;
            -k1*x(1) -k2*x(2) -k3*x(3) k4*x(4) 0 0 0 k8*x(8) 0 0;
            0 0 0 -k4*x(4) k5*x(5) k6*x(6) 0 0 -k9*x(9) 0;
            0 0 k3*x(3) 0 0 0 k7*x(7) 0 k9*x(9) -k10*x(10)
        ];
        J(8,:) = J(8,:)/k1;
        J(9,:) = J(9,:)/k1;
        J(10,:) = J(10,:)/k1;
end

function [C] = dJ_dw()
    a = 2;
    g = 9.81;
    w = 8/(g*pi^2);
    k1 = 0.0215*w*300/(0.255^5)*a;
    k2 = 0.0230*w*150/(0.255^5)*a;
    k3 = 0.0287*w*150/(0.15^5)*a;
    k4 = 0.0237*w*150/(0.15^5)*a;
    k5 = 0.0234*w*300/(0.205^5)*a;
    k6 = 0.0239*w*150/(0.205^5)*a;
    k7 = 0.0229*w*150/(0.205^5)*a;
    k8 = 0.0215*w*150/(0.255^5)*a;
    k9 = 0.0263*w*300/(0.205^5)*a;
    k10 = 0.0234*w*300/(0.205^5)*a;
    C = zeros([10 10]);
    C(8,:) = [-k1 -k2 -k3 k4 0 0 0 k8 0 0]/k1;
    C(9,:) = [0 0 0 -k4 k5 k6 0 0 -k9 0]/k1;
    C(10,:) = [0 0 k3 0 0 0 k7 0 k9 -k10]/k1;
end 