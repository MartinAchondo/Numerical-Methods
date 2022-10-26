
!Programa que resuleve sistema de edo con Runge Kutta 4
program main
    implicit none
    integer,parameter :: n=2
    integer :: i,mm
    real(8) :: tf,tout,ti,dt,ai(n),alpha,beta
    real(8), allocatable, dimension(:) :: tp,tp2
    real(8), allocatable, dimension(:,:) :: ap

    ti = 0
    ai(1) = 0
    ai(2) = 1    

    !Se lee de archivo tf, alpha y beta
    call ReadFile(tf,alpha,beta)

    dt = 0.01
    tout = 0.03
    mm = floor(tf/tout) + 2

    allocate (tp(mm))
    allocate (ap(n,mm))
    allocate (tp2(mm))

    !Se corre RK4
    call Main_RK4(n,ti,tf,tout,dt,ai,ap,tp,mm,alpha,beta)

    !write(*,*) tp
    write(*,*)
    !write(*,*) ap(1,:)

    do i=1,mm
        tp2(i) = tp(i) + log(beta)
    end do

    !Se escribe en archivo de salida
    call WriteFile(mm,tp,ap(1,:),tp2)


end program main


! Edos a resolver
subroutine Derivs(x,y,dy,n,a,b)

    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: x,y(n),a,b
    real(8), intent(out) :: dy(n)
    real(8) :: k

    k = 1/(b*exp(x))

    dy(1) = y(2)
    dy(2) = -3*y(2) - (2-a*k*(1-k**2))*y(1)

end subroutine


!Main de RK4
subroutine Main_RK4(n,xi,xf,xout,dx,yi,yp,xp,mm,a,b)

    implicit none
    integer, intent(in) :: n,mm
    real(8), intent(in) :: xi,xf,dx,yi(n),xout,a,b
    real(8), intent(out) :: yp(n,mm),xp(mm)
    integer :: m,i
    real(8) :: x,xend,h,y(n)

    x = xi
    m = 1
    xp(m) = x

    do i=1,n
        yp(i,m) = y(i)
        y(i) = yi(i)
    end do

    do while(.true.)
        xend = x + xout
        if(xend.gt.xf) then
            xend = xf
        end if
        h = dx
        call Integrator(x,y,n,h,xend,a,b)
        m = m + 1
        xp(m) = x
        do i=1,n
            yp(i,m) = y(i)
        end do
        if(x.ge.xf) then
            exit
        end if
    end do

end subroutine


subroutine Integrator(x,y,n,h,xend,a,b)
    
    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: a,b
    real(8), intent(inout) :: x,y(n),xend,h

    do while(.true.)

        if((xend-x).lt.h) then
            h = xend - x
        end if
        call RK4(x,y,n,h,a,b)
        if(x.ge.xend) then
            exit
        end if

    end do

end subroutine


subroutine RK4(x,y,n,h,a,b)

    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: h,a,b
    real(8), intent(inout) :: x,y(n)
    real(8) :: ye(n),ym(n),sl(n),k1(n),k2(n),k3(n),k4(n)

    call Derivs(x,y,k1,n,a,b)
    ym(:) = y(:) + k1(:)*h/2
    
    call Derivs(x+h/2,ym,k2,n,a,b)
    ym(:) = y(:) + k2(:)*h/2

    call Derivs(x+h/2,ym,k3,n,a,b)
    ye(:) = y(:) + k3(:)*h

    call Derivs(x+h,ye,k4,n,a,b)
    sl(:) = (k1(:) + 2*(k2(:)+k3(:)) + k4(:))/6
    y(:) = y(:) + sl(:)*h

    x = x + h

end subroutine


! Para escribir archivo salida
subroutine WriteFile(n,x,y,z)
    implicit none
    integer, intent(in) :: n
    integer :: i,stat
    real(kind=8), intent(in) :: x(n),y(n),z(n) 

    open(unit=20,file='Results_data.txt',iostat=stat,action='write')

        if(stat/=0) then
            write(*,*)'Error al abrir el archivo con iostat',stat
        end  if

        do i=1,n
            write(20,'(3F12.6)')x(i),y(i),z(i)
        end do
    close(unit=20,iostat=stat)

    if(stat/=0) then
        write(*,*)'Error al cerrar el archivo con iostat',stat
    end if

    write(*,*)
    write(*,*) "Archivo Guardado"
    write(*,*)

end subroutine WriteFile


! Para leer archivo de entrada
subroutine ReadFile(tf,a,b)
    implicit none
    integer :: stat
    real(8), intent(inout) :: tf,a,b

    open(unit=10,file='input_data.txt',iostat=stat,action='read')

        if(stat/=0) then
            write(*,*)'Error al abrir el archivo con iostat',stat
        end  if

        read(10,*)tf
        read(10,*)a
        read(10,*)b

    close(unit=10,iostat=stat)

    if(stat/=0) then
        write(*,*)'Error al cerrar el archivo con iostat',stat
    end if

    write(*,*)
    write(*,*) "Archivo Leido:"
    write(*,*)
    write(*,*) 'Tiempo de Integracion', tf
    write(*,*) 'Alpha', a
    write(*,*) 'Beta', b
    write(*,*)

end subroutine ReadFile


