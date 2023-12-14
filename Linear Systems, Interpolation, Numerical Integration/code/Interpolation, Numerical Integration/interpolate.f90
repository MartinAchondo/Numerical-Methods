
program interpolate
    implicit none
    integer, parameter :: n=25
    real(8) :: r(n),y(n),a(n),b(n),c(n),d(n),r2(n-1),v(n-1),x(2*n-1),yf(2*n-1)
    real(8) :: eval_spline
    integer :: i,z

    call ReadFile(r,y,n,"data_ry.txt")

    a = y

    call spline(n,r,a,b,c,d)

    z = 1
    do i=1,n-1
        r2(i) = (r(i+1) + r(i))/2
        v(i) = eval_spline(n,a,b,c,d,r,r2(i))
        x(z) = r(i)
        x(z+1) = r2(i)
        yf(z) = y(i)
        yf(z+1) = v(i)
        z = z + 2
    end do
    x(2*n-1) = r(n)
    yf(2*n-1) = y(n)

    write(*,*)
    write(*,*) 'Valores Interpolados r,y  (tabla 2)'
    do i=1,2*n-1
        write(*,*) x(i),yf(i)
    end do

    call WriteFile(x,yf,2*n-1,"data_ry_inter.txt")

    


end program interpolate


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
    character(11),intent(in) :: file
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
            write(20,'(2F15.10)')a(i),b(i)
        end do
    close(unit=20,iostat=stat)

    if(stat/=0) then
        write(*,*)'Error al cerrar el archivo con iostat',stat
    end if

    write(*,*)
    write(*,*) "Archivo Guardado"
    write(*,*)

end subroutine WriteFile














