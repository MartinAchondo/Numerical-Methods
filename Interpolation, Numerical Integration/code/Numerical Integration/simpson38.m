function [val1]= simpson38(h,f1,f2,f3,f4);
    val1= 3*h*((f1+3*(f2+f3)+f4)/8);
end 
    