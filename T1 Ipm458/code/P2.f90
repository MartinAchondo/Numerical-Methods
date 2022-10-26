
!Programa para resolver un sistema no lineal de ecuaciones
!utilizando el metodo de newton y broyden

program P2
    implicit none
    real(kind=8) :: x0(2),xreal(2),xf(2),es
    real(kind=8) :: a

    es = epsilon(a)

    x0 = (/-0.5,1.4/)
    xreal = (/0d0,1d0/)

    write(*,*)
    write(*,*)
    write(*,*) '======= METODO DE NEWTON ======='
    write(*,*)

    call newton(es,x0,xreal,xf)
    write(*,*)
    write(*,*) 'Solucion:'
    write(*,'(3F16.10)') xf

    write(*,*)
    write(*,*)
    write(*,*) '======= METODO DE BROYDEN ======='
    write(*,*)

    call broyden(es,x0,xreal,xf)
    write(*,*)
    write(*,*) 'Solucion:'
    write(*,'(3F16.10)') xf

end program P2

! subrutina para evaluar la funcion vecotrial
subroutine f_k(x,f)
    real(kind=8)::x(2)
    real(kind=8) :: f(2)
    f(1) = (x(1)+3)*(x(2)**3-7)+18
    f(2) = sin(x(2)*exp(x(1))-1)
end subroutine

!subrutina para evaluar la matriz jacobiana
subroutine jacb(x,J)
    real(kind=8) ::x(2)
    real(kind=8) :: J(2,2)
    J(1,1) = x(2)**3-7
    J(1,2) = (3*x(2)**2)*(x(1)+3)
    J(2,1) = cos(x(2)*exp(x(1))-1)*exp(x(1))*x(2)
    J(2,2) = cos(x(2)*exp(x(1))-1)*exp(x(1))
end subroutine

!subrutina para calcular el error norma2
subroutine norma2(x,xreal,n,norm)
    integer :: n,i
    real(kind=8) :: norm,n1,n2
    real(kind=8) :: x(n),xreal(n) 
    n1 = 0
    n2 = 0
    do i=1,n
        n1 = n1 + (x(i)-xreal(i))**2
        n2 = n2 + xreal(i)**2
    end do
    norm = sqrt(n1/n2)
end subroutine

!subrutina para calcualr la norma de un vector
subroutine norma_vec(x,n,norm)
    integer :: n,i
    real(kind=8) :: norm,n1
    real(kind=8) :: x(n)
    n1 = 0
    do i=1,n
        n1 = n1 + x(i)**2
    end do
    norm = n1
end subroutine

!subrutina para utilizar el metodo de newton
subroutine newton(es,x0,xreal,xf)
    implicit none
    integer(kind=4), parameter :: imax=500,n=2
    real(kind=8)::x0(n),xreal(n),xf(n)
    integer(kind=4) :: iter
    real(kind=8) :: xr(n),xrold(n),es,b(n),Acopy(n,n),norm,fr(n),J(n,n),normv,norm2
    xrold = x0
    xr = xrold
    iter = 0
    write(*,*) "      Iteracion , Vector Solucion :  x1 ,  x2  , Error norma2 real  , Error norma2 aproximado"
    do while (iter.lt.imax)
        xrold = xr
        iter = iter + 1
        call f_k(xrold,fr)
        b = -fr
        call jacb(xrold,J)
        Acopy = J
        call sist_lineal_2d(Acopy,b)
        xr = xrold + b
        call norma2(xr,xreal,n,norm)
        call norma_vec(xr,n,normv)
        if(normv.ne.0) then
            call norma2(xrold,xr,n,norm2)
        end if
        write(*,*) iter,xr(1),xr(2),norm,norm2
        if(norm.lt.es.or.iter.gt.imax) then
            exit
        end if
    end do
    xf = xr
end subroutine

!subrutina para utilizar el metodo de broyden
subroutine broyden(es,x0,xreal,xf)
    implicit none
    integer(kind=4), parameter :: imax=500,n=2
    real(kind=8)::x0(n),xreal(n),xf(n)
    integer(kind=4) :: iter,i,j
    real(kind=8) :: xr(n),xrold(n),es,b(n),Acopy(n,n),norm,fr(n),yk(n),Bk(n,n),frk(n),sk(n),fksk(n,n),skk,normv,norm2
    xrold = x0
    xr = xrold
    iter = 0
    write(*,*) "      Iteracion , Vector Solucion :  x1 ,  x2  , Error norma2 real  , Error norma2 aproximado"
    write(*,*)
    call jacb(xrold,Bk)

    do while (iter.lt.imax)
        xrold = xr
        iter = iter + 1
        call f_k(xrold,fr)
        b = -fr
        sk = b
        Acopy = Bk
        call sist_lineal_2d(Acopy,sk)
        xr = xrold + sk
        call norma2(xr,xreal,n,norm)
        call norma_vec(xr,n,normv)
        if(normv.ne.0) then
            call norma2(xrold,xr,n,norm2)
        end if
        write(*,*) iter,xr(1),xr(2),norm,norm2
        if(norm.lt.es.or.iter.gt.imax) then
            exit
        end if

        call f_k(xr,frk)

        yk = frk - fr
        skk = 0
        do i=1,n
            do j=1,n
                fksk(i,j) = (frk(i))*sk(j)
            end do
            skk = skk + sk(i)**2
        end do

        Bk = Bk +(fksk)/skk

    end do
    xf = xr

end subroutine

! subrutina pararesolver el sistema lineal utilizando lapack
subroutine sist_lineal_2d(A,b)
    real(kind=8) :: A(2,2),b(2)
    integer :: IPIV(2),NRHS,LDA,LDB,INFO
    NRHS = 1 
    LDA = 2 
    LDB = 2 
    call DGESV(2,NRHS,A,LDA,IPIV,b,LDB,INFO)
end subroutine