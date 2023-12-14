
program inverse_inter
    implicit none
    integer,parameter :: n=49
    real(8) :: r(n),y(n),AA,RR,rn,h,e,F,M
    real(8) :: a(n),b(n),c(n),d(n)
    integer :: i
    real(8) :: sigma(n),rx,eval_spline

    call ReadFile(r,y,n,"data_ry_inter.txt")
    AA = 1189.0921636889593
    RR = 46.380989762466697
    rn = 42.338586884955298
    F = 20000
    h = 6
    e = RR - rn
    M = F*(RR + h)

    do i=1,n
        sigma(i) = F/AA + M*(rn-r(i))/(AA*e*r(i))
        write(*,*) r(i),sigma(i)
    end do

    call WriteFile(r,sigma,n,"data_rs_inter.txt")

    a = sigma
    call spline(n,r,a,b,c,d)
    
    call findroot(n,r,a,b,c,d,rx)

    write(*,*) 'Valor de: r , sigma(r)'
    write(*,*)rx,eval_spline(n,a,b,c,d,r,rx)
    write(*,*)
    
end program inverse_inter

subroutine findroot(n,r,a,b,c,d,xf)
    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: r(n),a(n),b(n),c(n),d(n)
    real(8), intent(out) :: xf
    integer :: i,index
    real(8) :: v1,K(4),w,tol

    v1 = a(1)
    do i=2,n
        if (v1*a(i).lt.0) then
            index = i - 1
            exit
        end if
    end do

    K = (/a(index),b(index),c(index),d(index)/)
    w = r(index)

    tol = 1.0d-9
    call newton_rapshon(K,w,4,w,tol,xf)

end subroutine findroot


subroutine newton_rapshon(K,w,n,x0,tol,xf)
    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: K(n),w,x0,tol
    real(8), intent(out) :: xf
    integer :: iter
    integer, parameter :: imax = 500
    real(8) :: xrold,xr,ea,fr,eval_f_df
    xrold = x0
    xr = xrold
    iter = 0
    do while (iter.lt.imax)
        xrold = xr
        iter = iter + 1
        fr = eval_f_df(xrold,K,w,n)
        xr = xrold - fr
        if (xr.ne.0) then
            ea = abs((xr-xrold)/xr)*100
        end if
        if(ea.lt.tol.or.iter.gt.imax) then
            exit
        end if
    end do
    write(*,*) 'Newton Rapshon resuelto: Iteracion, ea'
    write(*,*) iter,ea
    write(*,*)
    xf = xr
end subroutine

real(8) function eval_f_df(x,K,w,n)
    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: K(n),x,w
    integer :: j
    real(8)::sum1,sum2
    sum1 = 0
    do j = 1,n
        sum1 = sum1 + K(j)*(x-w)**(j-1)
    end do

    sum2 = 0
    do j = 1,n-1
        sum2 = sum2 + (j)*K(j+1)*(x-w)**(j-1)
    end do
    
    eval_f_df = sum1/sum2

end function


real(8) function eval_spline(n,a,b,c,d,x,v)
    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: a(n),b(n),c(n),d(n),x(n),v
    integer :: i
    real(8) :: s,z

    do i=1,n
        if (x(i).le.v) then
            if (x(i+1).ge.v) then
                z = v - x(i)
                s = a(i) + b(i)*z + c(i)*z**2 + d(i)*z**3
            end if 
        end if
    end do
    eval_spline = s

end function eval_spline

subroutine spline(n,x,a,b,c,d)
    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: x(n),a(n)
    real(8), intent(out) :: b(n),c(n),d(n)
    integer :: i,j
    real(8) :: h(n),alfa(n),l(n),mu(n),z(n)

    do i=1,n-1
        h(i) = x(i+1) - x(i)
    end do

    do i=2,n-1
        alfa(i) = (3d0/h(i))*(a(i+1)-a(i))-((3d0/h(i-1))*(a(i)-a(i-1)))
    end do

    l(1) = 1d0
    mu(1) = 0d0
    z(1) = 0d0

    do i=2,n-1
        l(i) = 2d0*(x(i+1)-x(i-1))-h(i-1)*mu(i-1)
        mu(i) = h(i)/l(i)
        z(i) = (alfa(i)-h(i-1)*z(i-1))/l(i)
    end do

    l(n) = 1d0
    z(n) = 0d0
    c(n) = 0d0

    do j=n-1,1,-1
        c(j) = z(j)-mu(j)*c(j+1)
        b(j) = (a(j+1)-a(j))/h(j) - h(j)*(c(j+1)+2d0*c(j))/3d0
        d(j) = (c(j+1)-c(j))/(3d0*h(j))
    end do

end subroutine spline

subroutine ReadFile(a,b,n,file)
    implicit none
    character(17),intent(in) :: file
    integer, intent(in) :: n
    integer :: i,stat
    real(kind=8) :: a(n),b(n)

    open(unit=10,file=file,iostat=stat,action='read')

        if(stat/=0) then
            write(*,*)'Error al abrir el archivo con iostat',stat
        end  if

        do i=1,n
            read(10,*)a(i),b(i)
        end do
    close(unit=10,iostat=stat)

    if(stat/=0) then
        write(*,*)'Error al cerrar el archivo con iostat',stat
    end if

    write(*,*)
    write(*,*) "Archivo Leido:"
    do i=1,n
        write(*,*)a(i),b(i)
    end do
    write(*,*)

end subroutine ReadFile

subroutine WriteFile(a,b,n,file)
    implicit none
    character(17),intent(in) :: file
    integer, intent(in) :: n
    integer :: i,stat
    real(kind=8) :: a(n),b(n)

    open(unit=20,file=file,iostat=stat,action='write')

        if(stat/=0) then
            write(*,*)'Error al abrir el archivo con iostat',stat
        end  if

        do i=1,n
            write(20,'(2F10.4)')a(i),b(i)
        end do
    close(unit=20,iostat=stat)

    if(stat/=0) then
        write(*,*)'Error al cerrar el archivo con iostat',stat
    end if

    write(*,*)
    write(*,*) "Archivo Guardado"
    write(*,*)

end subroutine WriteFile