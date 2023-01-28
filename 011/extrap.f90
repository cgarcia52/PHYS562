program extrap_test

    use numtype 
    use thiele_approx
    implicit none
    real(dp), dimension(maxpt) :: zn, fn, an
    integer :: n, i
    real(dp), external :: func
    real(dp) :: x, dx, yy

    n = 25
    dx = 0.5_dp
    do i = 1, n
        zn(i) = i*dx
        fn(i) = func(zn(i))
        print *, zn(i),fn(i)
    end do
    print *, '--------'

    call thiele_coef(n,zn,fn,an)

    x = -1._dp + 1.e-30_dp
    yy = thiele_cf(x,n,zn,an)
    print *, x, func(x), yy

end program extrap_test

function func(x) result(f)

    use numtype, only : dp
    implicit none
    real(dp) :: x, f

    !f = x*cos(x) ! f = 1/(1+x**2) - x**3

    f = sin(x)/x

end function func


module thiele_approx

    use numtype
    implicit none
    integer, parameter :: maxpt = 50

    contains 

        subroutine theiel_coef(nn,zn,fn,an)

            use numtype
            implicit none
            real(dp), dimension(maxpt) :: zn, fn, an
            real(dp), dimension(maxpt,maxpt) :: gn
            integer :: nn, n, nz

            gn(1,:nn) = fn(1:nn)
            do n = 2, nn
                do nz = n, nn
                    gn(n,nz) = (gn(n-1, n-1)- gn(n-1,nz)) / & ((zn(nz)-zn(n-1)) * gn(n-1,nz))
                end do
            end do
            forall(n=1:nn) an(n) = gn(n,n)

        end subroutine theiel_coef

        function theiel_coef (z,nn,zn,an) result(cfrac)

            use numtype
            implicit none
            real(dp) :: z
            real(dp), dimension(maxpt) :: zn, an
            integer :: nn, n
            real(dp) :: cf0(2), cf1(2), cf(2), cfrac

            cf0(1) = 0._dp; cf0(2) = 1._dp
            cf1(1) = an(1); cf1(2) = 1._dp
            do n =1, nn-1
                cf = cf1 + (z-zn(n)) *an(n+1) * cf0
                cf0 = cf1; cf1 = cf
            end do
            cfrac = cf(1)/cf(2)

        end function
end module thiele_approx
