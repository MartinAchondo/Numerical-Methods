
program integrate
    implicit none
    integer,parameter :: n=49
    real(8) :: x(n),y(n),fi(n),A,R,rn

    call ReadFile(x,y,n,"data_ry_inter.txt")

    write(*,*)
    write(*,*) 'Valores utilizando Trapsun (A,R,rn)'
    fi(:) = 2*y(:)
    call Trapsun(x,fi,n,A)
    write(*,*) 'A:',A
    fi(:) = 2*x(:)*y(:)
    call Trapsun(x,fi,n,R)
    R = R/A
    write(*,*) 'R:',R
    fi(:) = 2*y(:)/x(:)
    call Trapsun(x,fi,n,rn)
    rn = A/rn
    write(*,*) 'rn',rn

    write(*,*)
    write(*,*) 'Valores utilizando Uneven (A,R,rn)'
    fi(:) = 2*y(:)
    call Uneven(x,fi,n,A)
    write(*,*) 'A:',A
    fi(:) = 2*x(:)*y(:)
    call Uneven(x,fi,n,R)
    R = R/A
    write(*,*) 'R:',R
    fi(:) = 2*y(:)/x(:)
    call Uneven(x,fi,n,rn)
    rn = A/rn
    write(*,*) 'rn',rn
    write(*,*)

end program integrate


subroutine Trapsun(x,y,n,Int)
        implicit none
        integer, intent(in) :: n
        real(8), intent(in) :: x(n),y(n)
        real(8), intent(out) :: Int
        integer :: i
        real(8) :: sum

        sum = 0
        do i=2,n
            sum = sum + (x(i)-x(i-1))*(y(i)+y(i-1))/2
        end do
        Int = sum

end subroutine Trapsun

subroutine Uneven(x,f,n,Int)
    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: x(n),f(n)
    real(8), intent(out) :: Int
    integer :: k,j
    real(8) :: h,hf,sum,Trap,Simp13,Simp38

    h = x(2) - x(1)
    k = 1
    sum = 0
    do j=2,n
        hf = x(j+1) - x(j)
        if (abs(h-hf).lt.0.000001) then
            if (k.eq.3) then
                sum = sum + Simp13(h,f(j-3),f(j-2),f(j-1))
                k = k - 1
            else
                k = k + 1
            end if
        else
            if (k.eq.1) then
                sum = sum + Trap(h,f(j-1),f(j))
            else
                if (k.eq.2) then
                    sum = sum + Simp13(h,f(j-2),f(j-1),f(j))
                else
                    sum = sum + Simp38(h,f(j-3),f(j-2),f(j-1),f(j))
                end if
                k = 1
            end if
        end if
        h = hf
    end do
    Int = sum

end subroutine Uneven


real(8) function Trap(h,f0,f1)
    real(8), intent(in) :: h,f0,f1
    Trap = h*(f0 + f1)/2
end function Trap

real(8) function Simp13(h,f0,f1,f2)
    real(8), intent(in) :: h,f0,f1,f2
    Simp13 = 2*h*(f0 + 4*f1 + f2)/6
end function Simp13

real(8) function Simp38(h,f0,f1,f2,f3)
    real(8), intent(in) :: h,f0,f1,f2,f3
    Simp38 = 3*h*(f0+3*(f1+f2)+f3)/8
end function Simp38



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