format compact
clear

A = readmatrix('problem2.xlsx','Range','A1:J10');
b = transpose(readmatrix('problem2.xlsx', 'Range', 'A13:J13'));





x0 = [-1;-1;-1.3;1;1.1;1;-1.2;-1;1;1]*0.4;
es = 0.0000001;

tic
[x,iter,r] = grad_conj(A,b,x0,es);
toc

tic
[x2,iter2,r2] = jacobi(A,b,x0,es);
toc

tic
[x3,iter3,r3] = gauss_seidel(A,b,x0,es);
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
