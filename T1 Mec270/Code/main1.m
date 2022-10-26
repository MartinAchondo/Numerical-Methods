format compact
clear

A = readmatrix('problem1.xlsx','Range','A1:N14');
b = transpose(readmatrix('problem1.xlsx', 'Range', 'A17:N17'));

x0 = [-1 ;-1 ;-1.3 ;1 ;1.1 ;1 ;-1.2 ;-1 ;1 ;1 ;1 ;1 ;1 ;1]*100;
es = 0.00001;

tic
[x,iter,r] = grad_conj(A'*A,A'*b,x0,es);
toc

tic
[x2,iter2,r2] = jacobi(A'*A,A'*b,x0,es);
toc

tic
[x3,iter3,r3] = gauss_seidel(A'*A,A'*b,x0,es);
toc

tic
[x4] = qr_method(A,b);
toc

tic
[x5] = cholesky_method(A'*A,A'*b);
toc

tic
[x6] = lu_method(A,b);
toc

