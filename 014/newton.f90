module setup

    use numtype
    implicit none
    integer, parameter :: nv = 3, lda = 5, nrhs = 1
    real(dp) :: a(nv) = (/10._dp, 28._dp, 8._dp/3   /)
    real(dp) :: x(lda), f(lda), deriv(lda, lda)

end module

program newton

    use setup
    implicit none

    real(dp) :: dx(lda, nrhs), ff
    integer :: info, i, ipiv(lda), maxstep
    real(dp), parameter :: eps = 1.e-10_dp

    x(1:nv) = (/40._dp, 60._dp, 50._dp /)

    maxstep = 15
    print '(7(5x,a))', 'x', ' y', ' z', '                       f_1', 'f_2', 'f_3', '|f|'

    do i = 1, maxstep

        call func(ff)
        print '(3f10.4,4x, 4e12.3)', x(1:nv),f(1:nv),ff
        if(ff <= eps) exit
        dx(1:nv,1) = -f(1:nv)

        info = 0
        call dgesv(nv,nrhs,deriv,lda,ipiv,dx,lda,info)
        if(info /= 0) stop ' info/= 0 '

        x(1:nv) = x(1:nv) + dx(1:nv, 1)

    end do
end program newton

subroutine func(ff)

    use setup
    implicit none

    real(dp) :: ff

    f(1) = a(1)*(x(2) - x(1))
    f(2) = a(2)*x(1) -x(2) - x(1)*x(3)
    f(3) = x(1)*x(2) - a(3)*x(3)
    ff = sqrt(dot_product(f(1:nv),f(1:nv)))

    deriv(1,1) = -a(1)
    deriv(2,1) = a(2) - x(3)
    deriv(3,1) = x(2)

    deriv(1,2) = a(1)
    deriv(2,3) = -1
    deriv(3,2) = x(1)

    deriv(1,3) = 0
    deriv(2,3) = -x(1)
    deriv(3,3) = a(3)

end subroutine func
