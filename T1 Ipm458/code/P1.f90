
!Programa para calcular las raices de un polinomio aplicado 
! a el estado de esfuerzo  para obtener los esfuerzos principales.

program P1
    implicit none
    real(kind=8)::C0,C1,C2
    real(kind=8) :: es,val(3),x0,K(4),kk(4),f,xx,kcopy(4),ea(3),eaf
    integer :: n,i,z

    es = 0.00005
    x0 = 1.1567

    call Constantes(C0,C1,C2)

    K = (/-C0,-C1,-C2,1d0/)  !coeficientes polinomio de partida
    n = size(K)

    call newton_rapshon(es,x0,K,n,0,xx,eaf) !primera raiz
    val(1) = xx
    ea(1) = eaf
    call divsin(K,n,val(1),kk)  !actualizacion de coeficientes
    z = 1
    do i = 2,n-1   !loop para otras raices y actualizacion de coeficientes
        call newton_rapshon(es,x0,kk,n,z,xx,eaf)
        val(i) = xx
        ea(i) = eaf
        kcopy = kk
        call divsin(kcopy,n,val(i),kk)
        z = z + 1
    end do

    write(*,*)
    write(*,*) 'Verificacion:'
    write(*,*) '  Raiz ecuacion,         Funcion evaluada,             Error relativo aproximado porcentual'
    do i = 1,n-1
        call evalf(val(i),K,n,0,f)
        write(*,*) val(i),f,ea(i)
    end do
     
    write(*,*)
    write(*,*) 

    contains
    !subrutina para evaluar las constantes
        subroutine Constantes(C0,C1,C2)
            implicit none
            real(kind=8)::ox,oy,oz,txy,tyx,txz,tzx,tzy,tyz,C0,C1,C2
            ox = -6
            oy = 18
            oz = -12
            txy = 9
            tyx = txy
            txz = -15
            tzx = txz
            tyz = 6
            tzy = 6
            C2 = ox + oy + oz
            C1 = txy**2 + tyz**2 + tzx**2 - ox*oy - oy*oz - oz*ox
            C0 = ox*oy*oz + 2*txy*tyz*tzx  - ox*tyz**2 - oy*tzx**2 - oz*txy**2
        end subroutine
end program P1

!subrutina para evaluar un polinomio dependiendo de los coeficientes asociados (input)
subroutine evalf(x,K,n,z,f)
    implicit none
    real(kind=8)::x,sum,f
    integer::n,j,z
    real(kind=8), intent(in)::K(n)
    sum = 0
    do j = 1,n-z
        sum = sum + K(j+z)*x**(j-1)
    end do
    f = sum
end subroutine

!subrutina para evaluar la derivada un polinomio dependiendo de los coeficientes asociados (input)
subroutine evaldf(x,K,n,z,df)
    implicit none
    real(kind=8)::x,sum,df
    integer::n,j,z
    real(kind=8), intent(in)::K(n)
    sum = 0
    do j = 1,n-1-z
        sum = sum + (j)*K(j+1+z)*x**(j-1)
    end do
    df = sum
end subroutine

!subrutina para division sintetica
subroutine divsin(K,n,r,kk)
    integer::n,i
    real(kind=8)::K(n)
    real(kind=8) :: r,b(n),kk(n)
    b(n) = K(n)
    i = n-1
    do while(i.ge.1)
        b(i) = K(i) + r*b(i+1)
        i = i-1
    end do
    do i=1,n
        kk(i) = b(i)
    end do
end subroutine

!subrutina para el metodo de newton rapshon
subroutine newton_rapshon(es,x0,K,n,z,xx,eaf)
    implicit none
    real(kind=8)::xx
    integer(kind=4) :: iter
    integer(kind=4), parameter :: imax=500
    real(kind=8) :: x0,xr,xrold,es,ea,fr,f,df,eaf
    integer::n,z
    real(kind=8), intent(in)::K(n)
    xrold = x0
    xr = xrold
    iter = 0
    do while (iter.lt.imax)
        xrold = xr
        iter = iter + 1
        call evalf(xrold,K,n,z,f)
        call evaldf(xrold,K,n,z,df)
        fr = f/df
        xr = xrold - fr
        if (xr.ne.0) then
            ea = abs((xr-xrold)/xr)*100
        end if
        if(ea.lt.es.or.iter.gt.imax) then
            exit
        end if
    end do
    xx = xr
    eaf = ea
end subroutine







