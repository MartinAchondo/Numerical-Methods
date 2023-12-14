
% Metodo del trapecio corregido 
function [valor5] = trapecio_corregido(n,x,f,g,k,h,a,b)
    valorr= trapecio(n,x,f,k,h);
    fa=eval(subs(g,k,a));
    fb=eval(subs(g,k,b));
    valor5= valorr + ((h^2)/12)*(fa-fb);
end