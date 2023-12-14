
format compact

sigma=[-65 5*sqrt(3) 0; 5*sqrt(3) -75 0; 0 0 -52.7];

disp("PI")
tic
[lambda,w] = pi(sigma)
toc
disp("------")

disp("IPI")
tic
[lambda,w] = ipi(sigma,@qr_method_linearsyst)
toc
disp("------")

disp("QR-Hessenberg")
tic
[vp,VP] = qr_hessenberg(sigma)
toc

disp("Eig")

[VP,vp] = eig(sigma)