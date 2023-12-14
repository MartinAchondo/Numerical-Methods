program main
    implicit none

    integer :: n,i,j
    real(8) :: tol
    real(8), allocatable, dimension(:) :: b,x,aij,aijT
    integer, allocatable, dimension(:) :: rows,cols,rowsT,colsT
    
    n = 100000
    tol = 1.0d-9 

    allocate(b(n))
    allocate(x(n))
    allocate(rows(n+1))
    allocate(cols(4*n-4))
    allocate(aij(4*n-4))
    allocate(rowsT(n+1))
    allocate(colsT(4*n-4))
    allocate(aijT(4*n-4))

    call create_A_b_sparse(n,aij,rows,cols,b)

    x = 0.0d0
    
    write(*,*)
    write(*,*) 'Metodo jacobi'
    write(*,*) 'Para n  =',n
    write(*,*)

    call jacobi_method(n,aij,rows,cols,b,x,tol)

    write(*,*)
    ! write(*,*) 'x:'
    ! do i=1,n
    !     write(*,*) x(i)
    ! end do
    ! write(*,*)

end program main


subroutine jacobi_method(n,aij,rows,cols,b,x,tol)
    implicit none
    integer, intent(in) :: n,rows(n+1),cols(4*n-4)
    real(8), intent(in) :: b(n),tol,aij(4*n-4)
    real(8), intent(inout) :: x(n)

    integer :: iter,imax, rowsUL(n+1),colsUL(4*n-4)
    real(8) :: x_old(n),vv(n),d,aijUL(4*n-4)
    real(8) :: er

    x_old(:) = x(:)
    call create_UL_sparse(n,aijUL,rowsUL,colsUL,d)

    er = 0.9d0
    imax = 10000

    do while (er.gt.tol)

        call matmul_sparse(n,rowsUL,colsUL,aijUL,x_old,vv)
  
        x(:) = (1.0d0/d)*(vv(:) + b(:))

        call matmul_sparse(n,rows,cols,aij,x,vv)
        er = norm2(b - vv)


        if (iter.gt.imax) then
            exit
        end if
        iter = iter + 1
        x_old(:) = x 

    end do

    write(*,*) 'Iterations', iter
    write(*,*) 'Residual',  er
    write(*,*) 'Max value', maxval(b-vv)


end subroutine jacobi_method


subroutine create_UL_sparse(n,aij,rows,cols,d)
    implicit none

    integer, intent(in) :: n
    real(8), intent(out) :: aij(4*(n-1)),d
    integer, intent(out) :: rows(n+1),cols(4*(n-1))
    integer :: i,nn,nnz,nr,j1,j2

    nn = n/2

    d =3.0d0

    nnz = 1
    rows(1) = nnz

    do i=1,n
        nr = 4
        if (i.eq.1 .or. i.eq.n .or. i.eq.nn .or. i.eq.nn+1) then
            nr = 3
        end if

        j1 = rows(i)
        j2 = j1 + nr -1

        if (i.eq.1) then
            aij(j1+1) = -1.0d0
            aij(j1+2) = 0.5d0
            cols(j1) = 1
            cols(j1+1) = 2
            cols(j1+2) = n
        else if (i.eq.n) then
            aij(j1) = 0.5d0
            aij(j1+1) = -1.0d0
            cols(j1) = 1
            cols(j1+1) = n-1
            cols(j1+2) = n
        else if (i.eq.nn .or.i.eq.nn+1) then
            aij(j1) = -1.0d0
            aij(j1+2) = -1.0d0
            cols(j1) = i-1
            cols(j1+1) = i
            cols(j1+2) = i+1
        else if (i.lt.nn) then
            aij(j1) = -1.0d0
            aij(j1+2) = -1.0d0
            aij(j1+3) = 0.5d0
            cols(j1) = i-1
            cols(j1+1) = i
            cols(j1+2) = i+1
            cols(j1+3) = n+1-i
        else if (i.gt.nn+1) then
            aij(j1) = 0.5d0
            aij(j1+1) = -1.0d0
            aij(j1+3) = -1.0d0
            cols(j1) = n+1-i
            cols(j1+1) = i-1
            cols(j1+2) = i
            cols(j1+3) = i+1
 
        end if

        rows(i+1) = j2 + 1
   
    end do

    aij(:) = -aij(:)

end subroutine create_UL_sparse

subroutine create_A_b_sparse(n,aij,rows,cols,b)
    implicit none

    integer, intent(in) :: n
    real(8), intent(out) :: b(n),aij(4*(n-1))
    integer, intent(out) :: rows(n+1),cols(4*(n-1))
    integer :: i,nn,nnz,nr,j1,j2

    nn = n/2

    b = 1.5d0
    b(1) = 2.0d0
    b(n) = 2.0d0
    b(nn) = 1.0d0
    b(nn+1) = 1.0d0

    nnz = 1
    rows(1) = nnz

    do i=1,n
        nr = 4
        if (i.eq.1 .or. i.eq.n .or. i.eq.nn .or. i.eq.nn+1) then
            nr = 3
        end if

        j1 = rows(i)
        j2 = j1 + nr -1

        if (i.eq.1) then
            aij(j1) = 3.0d0
            aij(j1+1) = -1.0d0
            aij(j1+2) = 0.5d0
            cols(j1) = 1
            cols(j1+1) = 2
            cols(j1+2) = n
        else if (i.eq.n) then
            aij(j1) = 0.5d0
            aij(j1+1) = -1.0d0
            aij(j1+2) = 3.0d0
            cols(j1) = 1
            cols(j1+1) = n-1
            cols(j1+2) = n
        else if (i.eq.nn .or.i.eq.nn+1) then
            aij(j1) = -1.0d0
            aij(j1+1) = 3.0d0
            aij(j1+2) = -1.0d0
            cols(j1) = i-1
            cols(j1+1) = i
            cols(j1+2) = i+1
        else if (i.lt.nn) then
            aij(j1) = -1.0d0
            aij(j1+1) = 3.0d0
            aij(j1+2) = -1.0d0
            aij(j1+3) = 0.5d0
            cols(j1) = i-1
            cols(j1+1) = i
            cols(j1+2) = i+1
            cols(j1+3) = n+1-i
        else if (i.gt.nn+1) then
            aij(j1) = 0.5d0
            aij(j1+1) = -1.0d0
            aij(j1+2) = 3.0d0
            aij(j1+3) = -1.0d0
            cols(j1) = n+1-i
            cols(j1+1) = i-1
            cols(j1+2) = i
            cols(j1+3) = i+1
        end if
        rows(i+1) = j2 + 1  
    end do
end subroutine create_A_b_sparse


subroutine matmul_sparse(n,rows,cols,aij,v,x)
    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: v(n),aij(4*(n-1))
    integer, intent(in) :: rows(n+1),cols(4*(n-1))
    real(8), intent(out) :: x(n)
    integer :: row,col,j1,j2,j
    real(8) :: a

    do row=1,n
        j1 = rows(row)
        j2 = rows(row+1)-1; 

        x(row) = 0
        do j=j1,j2
           col = cols(j)
           a = aij(j)
           x(row) = x(row) + a*v(col) 
        end do
     end do  

end subroutine matmul_sparse
