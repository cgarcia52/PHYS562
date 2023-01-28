module setup

    use numtype
    implicit none

    integer, parameter :: n_eq = 3
    real(dp), parameter :: hbar = 1._dp, hbar2 = hbar**2, mass = 1._dp
    real(dp), parameter :: omega = 1._dp, x0 = sqrt(hbar/(mass*omega))

    real(dp) :: energy, xmax, dstep

    real(dp), allocatable, dimension(:,:) :: wf

    integer :: imax

end module setup


program hold

    use setup
    use chebyshev
    implicit none

    real(dp) :: emin, emax, e0
    real(dp), external :: psi0
    integer :: nch, iz

    xmax = 15
    dstep = 0.001_dp

    imax = abs(nint(xmax/dstep))

    allocate(wf(-imax:imax, 2))

    emin = 0
    emax = 5
    nch = 35

    call chebyex(psi0, nch, cheb, emin, emax)
    call chebyzero(nch, cheb, emin, emax, z0, iz0)

    print *, z0(1:iz0)

    do iz = 1, iz0

        e0 = z0(iz)
        print *, iz, e0
        call wavef(iz,e0)

    end do

end program hold

subroutine wavef(iw, e)

    use setup
    implicit none
    real(dp), intent(in) :: e

    real(dp) :: x, psi(n_eq), parity
    integer :: iw, i

    energy  = e
    x = xmax

    psi(1) = exp(-x**2/(2*x0**2))
    psi(2) = - x/x0**2 * psi(1)
    psi(3) = 0

    do while (x > 0)
        call rk4step(x, -dstep, psi)
    end do

    if (abs(psi(1)) > abs(psi(2))) then
        parity = 1
    else
        parity = -1
    end if

    x = xmax

    psi(1) = exp(-x**2/(2*x0**2)) /sqrt(2*psi(3))
    psi(2) = -x/x0**2 * psi(1)
    psi(3) = 0
    i = imax+1

    do while (x > 0)
        i = i-1
        wf(i,1) = x; wf(-i,1) = -x
        wf(i,2) = psi(1); wf(-i,2) = parity * psi(1)
        call rk4step(x, -dstep, psi)
    end do

    print *, ' || ', 2*psi(3)

    do i = -imax/2, imax/2
        write(unit = 20+iw, fmt = '(2f15.5)') wf(i,1), wf(i,2)
    end do

end subroutine wavef


function psi0(e)

    use setup
    implicit none
    real(dp), intent(in) :: e
    real(dp) :: psi0
    real(dp) :: x, psi(n_eq)

    energy = e
    x = xmax

    psi(1) = exp(-x**2/(2*x0**2))
    psi(2) = -x/x0**2 * psi(1)
    psi(3) = 0

    do while (x > 0)
        call rk4step(x, -dstep, psi)
    end do

    psi0 = psi(1) * psi(2)

end function psi0


subroutine rk4step(x,h,y)

    use setup
    implicit none
    real(dp), intent(inout) :: x
    real(dp), intent(in) :: h
    real(dp), intent(inout), dimension(n_eq) :: y
    real(dp), dimension(n_eq) :: k1, k2, k3, k4, dy

    k1 = kv(x,      h, y)
    k2 = kv( x+h/2, h, y+k1/2)
    k3 = kv( x+h/2, h, y+k2/2)
    k4 = kv( x+h,   h, y+k3)

    dy = ( k1 + 2*k2 + 2*k3 + k4) / 6

    y = y + dy
    x = x + h

    contains

        function kv ( x, dx, y) result(k)

            use setup
            real(dp), intent(in) :: x
            real(dp), intent(in) :: dx
            real(dp), intent(in), dimension(n_eq) :: y
            real(dp), dimension(n_eq) :: f, k

            f(1) = y(2)
            f(2) = - 2*mass/hbar2 * ( energy - potential(x) ) * y(1)
            f(3) = -y(1)**2
            k = dx * f

        end function kv

        function potential (x)

            use setup
            real(dp) :: x, potential

            potential = 0.5_dp * mass * omega**2 * x**2

        end function potential

end subroutine rk4step