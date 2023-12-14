
program p1
    implicit none
    character(10) :: file1
    integer,parameter :: n = 16
    real(8),parameter :: tol = 1.0d-8
    integer :: er,o(n),q,i,oi(n)
    real(8) :: a(n,n),b(n),x(n),s(n),ai(n,n),a2(n,n),a3(n,n),e1,e2,e3,error_app,bi(n),az(n,n),si(n),ai1(n,n),xi(n),aco(n,n)
    real(8) :: a1(n,n),b1(n),x1(n),b2(n),x2(n),b3(n),x3(n),bc(n),bc2(n),bc3(n)

    file1 = "data_a.txt"
    write(*,*) file1
    call ReadFileA(a,n,file1)
    call ReadFileB(b,n,"dat_b1.txt")
    a1 = a
    a2 = a
    a3 = a
    aco = a
    ai1 = a
    b1 = b
    bc = b
    
    call ReadFileB(b,n,"dat_b2.txt")
    b2 = b
    bc2 = b2
    call ReadFileB(b,n,"dat_b3.txt")
    b3 = b
    bc3 = b3

    call Ludecomp(a1,b1,n,tol,x1,er,o,s)

    write(*,*)
    write(*,*) 'Resultado Caso sin Viento'
    do q=1,n
        write(*,*) x1(q)
    end do
    write(*,*)

    


    call Substitute(a1,o,n,b2,x2)

    write(*,*)
    write(*,*) 'Resultado Viento +x'
    do q=1,n
        write(*,*) x2(q)
    end do
    write(*,*)



    call Substitute(a1,o,n,b3,x3)

    write(*,*)
    write(*,*) 'Resultado Viento -x'
    do q=1,n
        write(*,*) x3(q)
    end do
    write(*,*)

    write(*,*)
    write(*,*) 'Resultado Lapcack sin viento'
    call sist_lineal(a,bc,n)

   do q=1,n
    write(*,*) bc(q)
   end do


   write(*,*)
   write(*,*) 'Resultado Lapcack viento +x'
   call sist_lineal(a2,bc2,n)

  do q=1,n
   write(*,*) bc2(q)
  end do


  write(*,*)
  write(*,*) 'Resultado Lapcack viento -x'
  call sist_lineal(a3,bc3,n)

 do q=1,n
  write(*,*) bc3(q)
 end do

 write(*,*)
 write(*,*)'Error Sin Viento:'
 write(*,*) error_app(x1,bc,n,e1)
 write(*,*) 'Error Viento +x'
 write(*,*) error_app(x2,bc2,n,e2)
 write(*,*) 'Error Viento -x'
 write(*,*) error_app(x3,bc3,n,e3)


 call InverseMa(ai1,n,tol,oi,si,er,bi,xi,ai)
 write(*,*)
 write(*,*) 'Matriz Inversa'
 do i=1,n
    write(*,'(16F8.4)')(ai(i,q),q=1,n)
 end do

 write(*,*)
 write(*,*) 'Comprobacion'
 az = matmul(aco,ai)
 do i=1,n
    write(*,'(16F8.4)')(az(i,q),q=1,n)
 end do

end program p1

real(8) function error_app(x1,x2,n,e)
    implicit none
    integer :: n,i,j
    real(8),intent(in) :: x1(n),x2(n)
    real(8) :: e,s1,s2
    s1 = 0
    s2 = 0
    do i=1,n
        s1 = s1 + (x1(i)-x2(i))**2
        s2 = s2 + x2(i)**2
    end do

    e = sqrt(s1/s2)

    error_app = e

end function error_app


subroutine Ludecomp(a,b,n,tol,x,er,o,s)
    implicit none
    integer :: n
    integer :: o(n),er
    real(8) :: a(n,n),b(n),x(n)
    real(8) :: s(n)
    real(8) :: tol
    er = 0
    call Decompose(a,n,tol,o,s,er)
    if (er.ne.-1) then
        call Substitute(a,o,n,b,x)
    end if

end subroutine Ludecomp

! Calcula inversa (LU debe estar hecha)
subroutine InverseM(a,n,o,er,b,x,ai)
    implicit none
    integer :: o(n),er,i,j,n
    real(8) :: a(n,n),b(n),x(n),ai(n,n)

    er = 0
    if (er.eq.0) then
        do i=1,n
            do j=1,n
                if (i.eq.j) then
                    b(j) = 1
                else
                    b(j) = 0
                end if
            end do
            call Substitute(a,o,n,b,x)
            do j=1,n
                ai(j,i) = x(j)
            end do
        end do
    end if
    
end subroutine InverseM

subroutine InverseMa(a,n,tol,o,s,er,b,x,ai)
    implicit none
    integer :: n
    integer :: o(n),er,i,j
    real(8) :: a(n,n),b(n),x(n),ai(n,n)
    real(8) :: s(n),tol

    er = 0
    call Decompose(a,n,tol,o,s,er)
    if (er.eq.0) then
        do i=1,n
            do j=1,n
                if (i.eq.j) then
                    b(j) = 1
                else
                    b(j) = 0
                end if
            end do
            call Substitute(a,o,n,b,x)
            do j=1,n
                ai(j,i) = x(j)
            end do
        end do
    end if
    
end subroutine InverseMa



subroutine Decompose(a,n,tol,o,s,er)
    implicit none
    integer :: i,j,k,o(n),er,n
    real(8) :: a(n,n),tol
    real(8) :: s(n),factor

    do i=1,n
        o(i) = i
        s(i) = abs(a(i,1))
        
        do j=2,n
            if (abs(a(i,j)).gt.s(i)) then
                s(i) = abs(a(i,j))
            end if            
        end do
    end do

    do k=1,n-1
        call Pivot(a,o,s,n,k)
        if (abs(a(o(k),k)/s(o(k))).lt.tol) then
            er = -1
            exit
        end if
        do i=k+1,n
            factor = a(o(i),k)/a(o(k),k)
            a(o(i),k) = factor
            do j=k+1,n
                a(o(i),j) = a(o(i),j) - factor*a(o(k),j)
            end do
        end do
    end do
    
    if (abs(a(o(k),k)/s(o(k))).lt.tol) then
        er = -1
    end if

end subroutine Decompose   


subroutine Pivot(a,o,s,n,k)
    implicit none
    integer :: ii,k,o(n),p,big,n
    real(8) :: a(n,n),dummy
    real(8) :: s(n)

    p = k
    big = abs(a(o(k),k)/s(o(k)))

    do ii=k+1,n
        dummy = abs(a(o(ii),k)/s(o(ii)))
        if (dummy.gt.big) then
            big = dummy
            p = ii
        end if
    end do

    dummy = o(p)
    o(p) = o(k)
    o(k) = dummy

end subroutine Pivot


subroutine Substitute(a,o,n,b,x)
    implicit none
    integer :: o(n),i,j,n
    real(8) :: a(n,n),b(n),x(n),sum

    do i=2,n
        sum = b(o(i))
        do j=1,i-1
            sum = sum - a(o(i),j)*b(o(j))
        end do
        b(o(i)) = sum
    end do

    x(n) = b(o(n))/a(o(n),n)

    do i=n-1,1,-1
        sum = 0
        do j=i+1,n
            sum = sum + a(o(i),j)*x(j)
        end do
        x(i) = (b(o(i))-sum)/a(o(i),i)
    end do

end subroutine Substitute


subroutine sist_lineal(A,b,n)
    implicit none
    integer, intent(in) :: n
    real(kind=8) :: A(n,n),b(n)
    integer :: IPIV(n),NRHS,LDA,LDB,INFO
    NRHS = 1 
    LDA = n 
    LDB = n 
    call DGESV(n,NRHS,A,LDA,IPIV,b,LDB,INFO)
end subroutine



subroutine ReadFileA(a,n,file)
    implicit none
    character(10) :: file
    integer :: i,j,stat,n
    real(8) :: a(n,n)

    open(unit=10,file=file,iostat=stat,action='read')

        if(stat/=0) then
            write(*,*)'Error al abrir el archivo con iostat',stat
        end  if

        do i=1,n
            read(10,*)(a(i,j),j=1,n)
        end do

    close(unit=10,iostat=stat)

    if(stat/=0) then
        write(*,*)'Error al cerrar el archivo con iostat',stat
    end if

    write(*,*) 'Archivo Leido'
    do i=1,n
        write(*,'(16F8.4)')(a(i,j),j=1,n)
    end do

end subroutine ReadFileA


subroutine ReadFileB(b,n,file)
    implicit none
    character(10) :: file
    integer :: i,j,stat,n
    real(8) :: b(n)

    open(unit=10,file=file,iostat=stat,action='read')

        if(stat/=0) then
            write(*,*)'Error al abrir el archivo con iostat',stat
        end  if

        do i=1,n
            read(10,*)b(i)
        end do

    close(unit=10,iostat=stat)

    if(stat/=0) then
        write(*,*)'Error al cerrar el archivo con iostat',stat
    end if

    write(*,*) 'Archivo Leido'
    do i=1,n
        write(*,*)b(i)
    end do

end subroutine ReadFileB

